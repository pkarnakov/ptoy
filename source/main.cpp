#include <iostream>

#include <GL/glew.h>
#include <SDL.h>
#include <SDL_opengl.h>
#include "game.hpp"
#include "logger.h"

#include <atomic>
#include <chrono>
#include <fstream>
#include <sstream>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

// TODO: reorder header includes (std first)
// TODO: fix formatting

using namespace std::chrono;

using Vect = vect; // TODO rename to Vect everywhere

namespace wrap {
auto glGetAttribLocation = [](GLuint prog, std::string varname) {
  GLint loc = ::glGetAttribLocation(prog, varname.c_str());
  fassert(loc != -1, "attribute '" + varname + "' not found");
  return loc;
};
auto glGetUniformLocation = [](GLuint prog, std::string varname) {
  GLint loc = ::glGetUniformLocation(prog, varname.c_str());
  fassert(loc != -1, "uniform '" + varname + "' not found");
  return loc;
};
void glCheckError() {
  GLenum e = glGetError();
  if (e != GL_NO_ERROR) {
    const std::string s(reinterpret_cast<const char*>(gluErrorString(e)));
    throw std::runtime_error(s);
  }
}

void glBufferDataReuse(
    GLenum target, GLsizeiptr size, const void* data, GLenum usage,
    GLuint vbo) {
  glBindBuffer(target, vbo);
  GLint cursize = 0;
  glGetBufferParameteriv(target, GL_BUFFER_SIZE, &cursize);

  if (cursize < size) { // FIXME: possible factor sizeof(GLfloat)
    glBufferData(target, size, data, usage);
  } else {
    glBufferSubData(target, 0, size, data);
  }
}

} // namespace wrap

std::string ReadFile(std::string path) {
  std::ifstream fs(path);
  fassert(fs.good(), "can't open file + '" + path + "'");
  std::stringstream ss;
  ss << fs.rdbuf();
  return ss.str();
}

void PrintShaderLog(GLuint shader) {
  if (glIsShader(shader)) {
    int max_length = 2048;
    int length = 0;

    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &max_length);

    char* log = new char[max_length];

    glGetShaderInfoLog(shader, max_length, &length, log);
    if (length > 0) {
      printf("%s\n", log);
    }

    delete[] log;
  } else {
    printf("Name %d is not a shader\n", shader);
  }
}

void PrintProgramLog(GLuint program) {
  if (glIsProgram(program)) {
    int max_length = 2048;
    int length = 0;

    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &max_length);

    char* log = new char[max_length];

    glGetProgramInfoLog(program, max_length, &length, log);
    if (length > 0) {
      printf("%s\n", log);
    }

    delete[] log;
  } else {
    printf("Name %d is not a program\n", program);
  }
}

GLuint CompileShaderSource(std::string src, GLenum type, std::string loc = "") {
  GLuint shader = glCreateShader(type);
  const GLchar* srcc[] = {src.c_str()};
  glShaderSource(shader, 1, srcc, NULL);
  glCompileShader(shader);

  GLint compiled = GL_FALSE;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
  if (compiled != GL_TRUE) {
    printf("Unable to compile shader %d in '%s'\n", shader, loc.c_str());
    PrintShaderLog(shader);
    fassert(false);
  }
  CHECK_ERROR();
  return shader;
}

GLuint CompileShaderFile(std::string path, GLenum type) {
  const std::string src = ReadFile(path);
  return CompileShaderSource(src, type, path);
}

GLuint CreateProgram(
    std::string vert_path, std::string frag_path, std::string geom_path) {
  GLuint program = glCreateProgram();

  std::vector<GLuint> todelete;

  if (vert_path.length()) {
    GLuint shader = CompileShaderFile(vert_path, GL_VERTEX_SHADER);
    glAttachShader(program, shader);
    todelete.push_back(shader);
    CHECK_ERROR();
  }

  if (frag_path.length()) {
    GLuint shader = CompileShaderFile(frag_path, GL_FRAGMENT_SHADER);
    glAttachShader(program, shader);
    todelete.push_back(shader);
    CHECK_ERROR();
  }

  if (geom_path.length()) {
    GLuint shader = CompileShaderFile(geom_path, GL_GEOMETRY_SHADER);
    glAttachShader(program, shader);
    todelete.push_back(shader);
    CHECK_ERROR();
  }

  glLinkProgram(program);
  { // error check
    GLint success = GL_TRUE;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (success != GL_TRUE) {
      PrintProgramLog(program);
      printf("Error linking program %d!\n", program);
      fassert(false);
    }
  }
  CHECK_ERROR();

  for (auto shader : todelete) {
    glDeleteShader(shader);
  }

  return program;
}

milliseconds last_frame_time;
Scal last_frame_game_time;
milliseconds last_report_time;
std::atomic<Scal> next_game_time_target;

std::unique_ptr<game> G;

GLfloat width, height;

std::array<Vect, 2> GetDomain() {
  Vect A(-1, -1);
  Vect B(-1 + 2. * width / kInitWidth, -1 + 2. * height / kInitHeight);
  return {A, B};
}

void SetUniformDomainSize(GLuint program) {
  auto dom = GetDomain();
  glUniform2ui(glGetUniformLocation(program, "screenSize"), width, height);
  glUniform2f(glGetUniformLocation(program, "domain0"), dom[0].x, dom[0].y);
  glUniform2f(glGetUniformLocation(program, "domain1"), dom[1].x, dom[1].y);
  CHECK_ERROR();
}

std::atomic<bool> flag_display;
std::atomic<bool> quit;
std::atomic<bool> pause;

using std::cout;
using std::endl;

struct RenderTask {
  GLuint program = 0;
  std::function<void()> render;
};

std::vector<RenderTask> tasks;

GLuint gTextureFrame = 0;

int frame_number;

void display() {
  if (flag_display) return;
  flag_display = true;

  const Scal fps = 30.0;
  milliseconds current_time = duration_cast<milliseconds>(
      high_resolution_clock::now().time_since_epoch());
  milliseconds time_past_from_last_frame = current_time - last_frame_time;

  auto new_frame_time = duration_cast<milliseconds>(
      high_resolution_clock::now().time_since_epoch());
  auto new_frame_game_time = G->PS->GetTime();

  const Scal frame_real_duration_s =
      (new_frame_time - last_frame_time).count() / 1000.;

  if ((new_frame_time - last_report_time).count() > 1000.) {
    std::cout << "fps: " << 1. / frame_real_duration_s << ", game rate="
              << (new_frame_game_time - last_frame_game_time) /
                     frame_real_duration_s
              << ", particles="
              << G->PS->GetNumParticles()
              //<< ", max_per_cell="
              //<< G->PS->GetNumPerCell()
              << ", t=" << G->PS->GetTime()
              << ", steps=" << G->PS->GetNumSteps() << std::endl;
    last_report_time = new_frame_time;
    // G->PS->Blocks.print_status();
  }

  const Scal game_rate_target = 1.;
  if (!pause) {
    next_game_time_target = new_frame_game_time + game_rate_target / fps;
  }

  last_frame_time = new_frame_time;
  last_frame_game_time = new_frame_game_time;
  // cout<<"Frame "<<frame_number++<<endl;

  for (auto& task : tasks) {
    task.render();
    CHECK_ERROR();
  }

  {
    std::lock_guard<std::mutex> lg(G->PS->m_buffer_);
    // G->R->DrawAll();
    G->PS->SetRendererReadyForNext(true);
  }

  // glFlush();

  /*
  glGenTextures(1, &gTextureFrame);
  glBindTexture(GL_TEXTURE_RECTANGLE, gTextureFrame);
  glCopyTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA, 0, 0, width, height, 0);
  */

  flag_display = false;
}

void cycle() {
  // omp_set_dynamic(0);
  // omp_set_nested(0);
  // omp_set_num_threads(std::thread::hardware_concurrency());
#ifdef _OPENMP
  // omp_set_num_threads(1);
#endif

#pragma omp parallel
  {
#pragma omp master
    std::cout << "Computation started, "
#ifdef _OPENMP
              << "OpenMP with " << omp_get_num_threads() << " threads"
#else
              << "single-threaded"
#endif
              << std::endl;
  }

  while (!quit) {
    G->PS->step(next_game_time_target, pause);
    // std::this_thread::sleep_for(milliseconds(50));
    if (pause) {
      std::this_thread::sleep_for(milliseconds(100));
    }
  }
  std::cout << "Computation finished" << std::endl;
}

// TODO: detect motionless regions

Vect GetDomainMousePosition(int x, int y) {
  Vect c;
  Vect A(-1., -1.),
      B(-1 + 2. * width / kInitWidth, -1. + 2. * height / kInitHeight);
  c.x = A.x + (B.x - A.x) * (x / width);
  c.y = B.y + (A.y - B.y) * (y / height);
  return c;
}

Vect GetDomainMousePosition() {
  int x, y;
  SDL_GetMouseState(&x, &y);
  return GetDomainMousePosition(x, y);
}

int main() {
  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    SDL_Log("Unable to initialize SDL: %s\n", SDL_GetError());
    return 1;
  }

  width = kInitWidth;
  height = kInitHeight;

  last_frame_game_time = 0.;
  pause = false;

  G = std::unique_ptr<game>(new game(width, height));

  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

  SDL_Window* window = SDL_CreateWindow(
      "ptoy", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height,
      SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

  if (window == NULL) {
    SDL_Log("Unable to create window: %s\n", SDL_GetError());
    return 1;
  }
  SDL_GLContext glcontext = SDL_GL_CreateContext(window);
  SDL_GL_SetSwapInterval(1);
  // SDL_GL_SwapWindow(window);

  glewInit();

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glEnable(GL_PROGRAM_POINT_SIZE);

  { // particles
    GLuint program = CreateProgram(
        "../shaders/vert.glsl", "../shaders/frag.glsl", "../shaders/geom.glsl");

    GLuint vbo_point;
    glGenBuffers(1, &vbo_point);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_point);
    GLint attr_point = wrap::glGetAttribLocation(program, "point");
    glEnableVertexAttribArray(attr_point);

    GLuint vbo_color;
    glGenBuffers(1, &vbo_color);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
    GLint attr_color = wrap::glGetAttribLocation(program, "color");
    glEnableVertexAttribArray(attr_color);

    GLuint tex_circle = 0;
    glGenTextures(1, &tex_circle);
    glBindTexture(GL_TEXTURE_2D, tex_circle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    const int w = 64;
    const int h = 64;
    std::vector<float> v(w * h);
    for (int j = 0; j < h; ++j) {
      for (int i = 0; i < w; ++i) {
        using std::sin;
        float fi = i;
        float fj = j;
        v[j * w + i] =
            (sin(0.2 + (fj / h) * 10 - sqr(fi / w - 0.5) * 10) + 1) * 0.5;
      }
    }
    glTexImage2D(
        GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, v.data());

    auto render = [program, vbo_point, vbo_color, tex_circle, &ps = G->PS,
                   attr_point, attr_color]() {
      glUseProgram(program);

      auto& particles = ps->GetParticles();
      std::vector<GLfloat> buf(particles.size() * 2);
      std::vector<GLfloat> buf_color(particles.size());
      {
        std::lock_guard<std::mutex> lg(ps->m_buffer_);
        for (size_t i = 0; i < particles.size(); ++i) {
          auto& p = particles[i];
          buf[i * 2] = p.p.x;
          buf[i * 2 + 1] = p.p.y;
          Scal f = 0.5 + p.v.length() / 3.; // color intensity
          f = std::min<Scal>(std::max<Scal>(f, 0.), 1.);
          buf_color[i] = f;
        }
      }

      wrap::glBufferDataReuse(
          GL_ARRAY_BUFFER, buf.size() * sizeof(GLfloat), buf.data(),
          GL_DYNAMIC_DRAW, vbo_point);
      glVertexAttribPointer(
          attr_point, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), NULL);

      wrap::glBufferDataReuse(
          GL_ARRAY_BUFFER, buf_color.size() * sizeof(GLfloat), buf_color.data(),
          GL_DYNAMIC_DRAW, vbo_color);
      glVertexAttribPointer(
          attr_color, 1, GL_FLOAT, GL_FALSE, 1 * sizeof(GLfloat), NULL);

      SetUniformDomainSize(program);
      glUniform1ui(
          glGetUniformLocation(program, "pointSize"),
          kRadius * kInitHeight / 2);

      CHECK_ERROR();
      glDrawArrays(GL_POINTS, 0, particles.size());
      CHECK_ERROR();
      glUseProgram(0);
    };

    tasks.push_back({program, render});
  }

  { // lines (frame, bonds, portals)
    GLuint program = CreateProgram(
        "../shaders/lines.vert.glsl", "../shaders/lines.frag.glsl",
        "../shaders/lines.geom.glsl");

    GLuint vbo_point;
    glGenBuffers(1, &vbo_point);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_point);
    GLint attr_point = wrap::glGetAttribLocation(program, "point");
    glEnableVertexAttribArray(attr_point);

    GLuint vbo_color;
    glGenBuffers(1, &vbo_color);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
    GLint attr_color = glGetAttribLocation(program, "color");
    glEnableVertexAttribArray(attr_color);

    GLuint vbo_width;
    glGenBuffers(1, &vbo_width);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_width);
    GLint attr_width = glGetAttribLocation(program, "width");
    glEnableVertexAttribArray(attr_width);

    auto render = [program, vbo_point, vbo_color, vbo_width, &ps = G->PS,
                   attr_point, attr_color, attr_width]() {
      glUseProgram(program);

      std::vector<GLfloat> buf;
      std::vector<GLfloat> buf_color;
      std::vector<GLfloat> buf_width;
      size_t nprim = 0;

      auto add = [&](Vect a, Vect b, Scal color, Scal width) {
        buf.push_back(a.x);
        buf.push_back(a.y);
        buf.push_back(b.x);
        buf.push_back(b.y);
        buf_color.push_back(color);
        buf_color.push_back(color);
        buf_width.push_back(width);
        buf_width.push_back(width);
        ++nprim;
      };
      auto rgb = [](float r, float g, float b) { return r; };

      { // draw portals
        const auto& portals = ps->GetPortals();
        const Scal kPortalWidth = 0.01;
        for (auto& pair : portals) {
          add(pair[0].begin, pair[0].end, rgb(0., .5, 1.), kPortalWidth);
          add(pair[1].begin, pair[1].end, rgb(1., .5, 0.), kPortalWidth);
        }
        if (ps->portal_stage_ == 0) {
          if (ps->portal_mouse_moving_) {
            add(ps->portal_begin_, ps->portal_current_, rgb(0., .5, 1.),
                kPortalWidth);
          }
        } else {
          add(ps->portal_prev_.first, ps->portal_prev_.second, rgb(0., .5, 1.),
              kPortalWidth);
          if (ps->portal_mouse_moving_) {
            add(ps->portal_begin_, ps->portal_current_, rgb(1., .5, 0.),
                kPortalWidth);
          }
        }
      }

      { // draw bonds
        const Scal kBondWidth = 0.02;
        const auto& nr = ps->GetNoRendering();
        const auto& pos = ps->GetBlockData().position;
        const auto& bbi = ps->GetBlockById();
        for (auto bond : ps->GetBonds()) {
          const auto& a = bbi[bond.first];
          const auto& b = bbi[bond.second];
          if (!nr.count(bond)) {
            add(pos[a.first][a.second], pos[b.first][b.second], 1, kBondWidth);
          }
        }
      }

      { // draw frozen
        // TODO: draw with texture or color of particles instead
        const Scal kFrozenWidth = kRadius * 2;
        const auto& pos = ps->GetBlockData().position;
        const auto& bbi = ps->GetBlockById();
        for (auto id : ps->GetFrozen()) {
          const auto& a = bbi[id];
          Vect c = pos[a.first][a.second];
          Vect d(kFrozenWidth * 0.25, 0);
          add(c - d, c + d, rgb(0., 0., 0.), kFrozenWidth);
        }
      }

      { // draw frame
        const Scal kFrameWidth = 0.005;
        auto dom = ps->GetDomain();
        Vect dom0 = dom.A;
        Vect dom1 = dom.B;
        add(Vect(dom0.x, dom0.y), Vect(dom1.x, dom0.y), 1, kFrameWidth);
        add(Vect(dom1.x, dom0.y), Vect(dom1.x, dom1.y), 1, kFrameWidth);
        add(Vect(dom1.x, dom1.y), Vect(dom0.x, dom1.y), 1, kFrameWidth);
        add(Vect(dom0.x, dom1.y), Vect(dom0.x, dom0.y), 1, kFrameWidth);
      }


      fassert_equal(buf.size(), nprim * 4);
      fassert_equal(buf_width.size(), nprim * 2);
      fassert_equal(buf_color.size(), nprim * 2);

      wrap::glBufferDataReuse(
          GL_ARRAY_BUFFER, buf.size() * sizeof(GLfloat), buf.data(),
          GL_DYNAMIC_DRAW, vbo_point);
      glVertexAttribPointer(
          attr_point, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), NULL);

      wrap::glBufferDataReuse(
          GL_ARRAY_BUFFER, buf_color.size() * sizeof(GLfloat), buf_color.data(),
          GL_DYNAMIC_DRAW, vbo_color);
      glVertexAttribPointer(
          attr_color, 1, GL_FLOAT, GL_FALSE, 1 * sizeof(GLfloat), NULL);

      wrap::glBufferDataReuse(
          GL_ARRAY_BUFFER, buf_width.size() * sizeof(GLfloat), buf_width.data(),
          GL_DYNAMIC_DRAW, vbo_width);
      glVertexAttribPointer(
          attr_width, 1, GL_FLOAT, GL_FALSE, 1 * sizeof(GLfloat), NULL);

      SetUniformDomainSize(program);
      CHECK_ERROR();
      glDrawArrays(GL_LINES, 0, nprim * 2);
      CHECK_ERROR();
      glUseProgram(0);
    };

    tasks.push_back({program, render});
  }

  std::thread computation_thread(cycle);

  G->PS->SetForce(Vect(0., 0.), false);

  // Main loop flag
  quit = false;

  // Event handler
  SDL_Event e;

  enum class MouseState { None, Force, Bonds, Pick, Freeze, Portal };
  MouseState mouse_state = MouseState::Force;

  size_t quit_count = 0;

  // While application is running
  while (!quit) {
    while (SDL_PollEvent(&e) != 0) {
      if (e.type == SDL_QUIT) {
        quit = true;
      } else if (e.type == SDL_KEYDOWN) {
        if (e.key.keysym.sym != 'q') {
          quit_count = 0;
        }
        switch (e.key.keysym.sym) {
          case 'q':
            if (++quit_count >= 3) {
              quit = true;
            } else {
              std::cout << "Quit counter: " << quit_count << "/3" << std::endl;
            }
            break;
          case 'n':
            mouse_state = MouseState::None;
            std::cout << "Mouse switched to None mode" << std::endl;
            break;
          case 'r':
            mouse_state = MouseState::Force;
            std::cout << "Mouse switched to Repulsive mode" << std::endl;
            G->PS->SetForceAttractive(false);
            break;
          case 'f':
            mouse_state = MouseState::Freeze;
            std::cout << "Mouse switched to Freeze mode" << std::endl;
            break;
          case 'p':
            mouse_state = MouseState::Pick;
            std::cout << "Mouse switched to Pick mode" << std::endl;
            break;
          case 'a':
            mouse_state = MouseState::Force;
            std::cout << "Mouse switched to Attractive mode" << std::endl;
            G->PS->SetForceAttractive(true);
            break;
          case 'b':
            mouse_state = MouseState::Bonds;
            std::cout << "Mouse switched to Bonds mode" << std::endl;
            break;
          case 'o':
            mouse_state = MouseState::Portal;
            std::cout << "Mouse switched to Portal mode" << std::endl;
            break;
          case 'i':
            std::cout << "Remove last portal" << std::endl;
            G->PS->RemoveLastPortal();
            break;
          case 'g':
            G->PS->InvertGravity();
            std::cout << (G->PS->GetGravity() ? "Gravity on" : "Gravity off")
                      << std::endl;
            break;
          case ' ':
            pause = !pause;
            break;
        }
      } else if (e.type == SDL_MOUSEMOTION) {
        switch (mouse_state) {
          case MouseState::Force:
            G->PS->SetForce(GetDomainMousePosition());
            break;
          case MouseState::Bonds:
            G->PS->BondsMove(GetDomainMousePosition());
            break;
          case MouseState::Pick:
            G->PS->PickMove(GetDomainMousePosition());
            break;
          case MouseState::Freeze:
            G->PS->FreezeMove(GetDomainMousePosition());
            break;
          case MouseState::Portal:
            G->PS->PortalMove(GetDomainMousePosition());
            break;
          case MouseState::None:
            break;
        }
      } else if (e.type == SDL_MOUSEBUTTONDOWN) {
        switch (mouse_state) {
          case MouseState::Force:
            G->PS->SetForce(GetDomainMousePosition(), true);
            break;
          case MouseState::Bonds:
            G->PS->BondsStart(GetDomainMousePosition());
            break;
          case MouseState::Pick:
            G->PS->PickStart(GetDomainMousePosition());
            break;
          case MouseState::Freeze:
            G->PS->FreezeStart(GetDomainMousePosition());
            break;
          case MouseState::Portal:
            G->PS->PortalStart(GetDomainMousePosition());
            break;
          case MouseState::None:
            break;
        }
      } else if (e.type == SDL_MOUSEBUTTONUP) {
        switch (mouse_state) {
          case MouseState::Force:
            G->PS->SetForce(false);
            break;
          case MouseState::Bonds:
            G->PS->BondsStop(GetDomainMousePosition());
            break;
          case MouseState::Pick:
            G->PS->PickStop(GetDomainMousePosition());
            break;
          case MouseState::Freeze:
            G->PS->FreezeStop(GetDomainMousePosition());
            break;
          case MouseState::Portal:
            G->PS->PortalStop(GetDomainMousePosition());
            break;
          case MouseState::None:
            break;
        }
      } else if (e.type == SDL_WINDOWEVENT) {
        switch (e.window.event) {
          case SDL_WINDOWEVENT_RESIZED:
            width = e.window.data1;
            height = e.window.data2;
            glViewport(0, 0, width, height);
            G->SetWindowSize(width, height);
            break;
        }
      }
    }

    auto gray = 0.5;
    glClearColor(gray, gray, gray, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    display();
    SDL_GL_SwapWindow(window);
  }

  // Finalize
  SDL_GL_DeleteContext(glcontext);
  SDL_DestroyWindow(window);
  SDL_Quit();

  computation_thread.join();

  return 0;
}
