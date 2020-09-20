#include <iostream>

#include <GL/glew.h>
#include <SDL.h>
#include <SDL_opengl.h>
#include "game.hpp"
#include "logger.h"

#include <atomic>
#include <chrono>
#include <fstream>
#include <map>
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
auto glGetAttribLocation = [](GLuint prog, std::string varname,
                              bool nofail = false) {
  GLint loc = ::glGetAttribLocation(prog, varname.c_str());
  if (!nofail) {
    fassert(loc != -1, "attribute '" + varname + "' not found");
  }
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

template <class T>
GLenum GetGlTypeEnum();

template <>
GLenum GetGlTypeEnum<GLfloat>() {
  return GL_FLOAT;
}
template <>
GLenum GetGlTypeEnum<GLint>() {
  return GL_INT;
}
template <>
GLenum GetGlTypeEnum<GLchar>() {
  return GL_BYTE;
}

template <class T, int ncomp>
struct VertexAttribute {
  VertexAttribute(std::string name, GLuint program, bool nofail=false) {
    glGenBuffers(1, &id);
    location = wrap::glGetAttribLocation(program, name, nofail);
    glEnableVertexAttribArray(location);
  }
  ~VertexAttribute() {
    glDeleteBuffers(1, &id);
  }
  void SetData(const std::vector<T>& data) const {
    wrap::glBufferDataReuse(
        GL_ARRAY_BUFFER, data.size() * sizeof(T), data.data(),
        GL_DYNAMIC_DRAW, id);
    glVertexAttribPointer(
        location, ncomp, GetGlTypeEnum<T>(), GL_FALSE, ncomp * sizeof(T),
        NULL);
  }
  void SetData(const std::vector<std::array<T, ncomp>>& data) const {
    wrap::glBufferDataReuse(
        GL_ARRAY_BUFFER, data.size() * sizeof(T) * ncomp, data.data(),
        GL_DYNAMIC_DRAW, id);
    glVertexAttribPointer(
        location, ncomp, GetGlTypeEnum<T>(), GL_FALSE, ncomp * sizeof(T),
        NULL);
  }
  void SetDataInt(const std::vector<std::array<T, ncomp>>& data) const {
    wrap::glBufferDataReuse(
        GL_ARRAY_BUFFER, data.size() * sizeof(T) * ncomp, data.data(),
        GL_DYNAMIC_DRAW, id);
    glVertexAttribIPointer(
        location, ncomp, GetGlTypeEnum<T>(), ncomp * sizeof(T), NULL);
  }
  GLuint id;
  GLint location;
};

struct Texture {
  Texture(GLenum target = GL_TEXTURE_2D) : target(target) {
    glGenTextures(1, &id);
    Bind();
    glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  }
  ~Texture() {
    glDeleteTextures(1, &id);
  }
  void SetData(
      const std::vector<float>& data, unsigned width, unsigned height) const {
    fassert_equal(data.size(), width * height);
    Bind();
    glTexImage2D(
        target, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, data.data());
  }
  void Bind() const {
    glBindTexture(target, id);
  }
  GLuint id;
  GLenum target;
};

struct Font {
  std::vector<float> data;
  size_t width;
  size_t height;
  float offset;
};

Font LoadFont(std::string path) {
  Font font;

  {
    const std::string p = path + ".size";
    std::ifstream f(p);
    fassert(f.good(), "can't open '" + p + "'");
    f >> font.width >> font.height >> font.offset;
  }

  {
    // file written by `font.py`
    const std::string p = path + ".data";
    std::ifstream f(p, std::ios::binary);
    fassert(f.good(), "can't open '" + p + "'");
    std::vector<unsigned char> buf(
        std::istreambuf_iterator<char>(f), {});
    fassert_equal(buf.size(), font.width * font.height);
    font.data.resize(font.width * font.height);
    for (size_t i = 0; i < font.data.size(); ++i) {
      font.data[i] = buf[i] / 255.;
    }
  }
  return font;
}

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

class Gui {
 public:
  using Color = std::array<float, 3>;

  static constexpr unsigned kButtonHeight = 50;
  static constexpr unsigned kButtonWidth = 100;

  struct Window {
    Vect size;
  };


  struct Button {
    std::string label;
    std::function<void()> handler;
    Color color;
    Vect lowcorner;
    Vect size;

    Button(std::string label, std::function<void()> handler, Color color)
        : label(label), handler(handler), color(color) {
      size = Vect(kButtonWidth, kButtonHeight);
    }
  };

  void AddButton(const Button& e) {
    buttons_.push_back(e);
  }
  const std::vector<Button>& GetButtons() const {
    return buttons_;
  }
  void SetWindowSize(const size_t width, const size_t height) {
    window_.size = Vect(width, height);
    UpdateButtonPositions();
  }
  void Update() {
    UpdateButtonPositions();
  }

 private:
  void UpdateButtonPositions() {
    Scal x = 0;
    for (auto& b : buttons_) {
      b.lowcorner = Vect(x, window_.size.y - kButtonHeight);
      x += kButtonWidth;
    }
  }

  std::vector<Button> buttons_;
  Window window_;
};

milliseconds last_frame_time;
Scal last_frame_game_time;
milliseconds last_report_time;
std::atomic<Scal> next_game_time_target;

std::unique_ptr<game> G;

Gui gui;

unsigned width, height;

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
std::vector<size_t> tasks_indices;

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

  for (size_t i : tasks_indices) {
    tasks[i].render();
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
  c.x = A.x + (B.x - A.x) * (float(x) / width);
  c.y = B.y + (A.y - B.y) * (float(y) / height);
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

  glClearColor(0.5, 0.5, 0.5, 1.0);
  glDrawBuffer(GL_FRONT);
  glClear(GL_COLOR_BUFFER_BIT);
  glDrawBuffer(GL_BACK);

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glEnable(GL_PROGRAM_POINT_SIZE);

  glProvokingVertex(GL_FIRST_VERTEX_CONVENTION);

  std::map<std::string, size_t> task2index;

  { // particles
    GLuint program = CreateProgram(
        "../shaders/vert.glsl", "../shaders/frag.glsl", "../shaders/geom.glsl");

    auto attr_point =
        std::make_shared<VertexAttribute<GLfloat, 2>>("point", program);
    auto attr_color =
        std::make_shared<VertexAttribute<GLfloat, 1>>("color", program);

    auto tex_circle = std::make_shared<Texture>();
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
    tex_circle->SetData(v, w, h);

    auto render = [program, tex_circle, &ps = G->PS, attr_point, attr_color]() {
      glUseProgram(program);

      auto& particles = ps->GetParticles();
      std::vector<std::array<GLfloat, 2>> buf(particles.size());
      std::vector<GLfloat> buf_color(particles.size());
      {
        std::lock_guard<std::mutex> lg(ps->m_buffer_);
        for (size_t i = 0; i < particles.size(); ++i) {
          auto& p = particles[i];
          buf[i] = {p.p.x, p.p.y};
          Scal f = 0.5 + p.v.length() / 3.; // color intensity
          f = std::min<Scal>(std::max<Scal>(f, 0.), 1.);
          buf_color[i] = f;
        }
      }

      tex_circle->Bind();

      attr_point->SetData(buf);
      attr_color->SetData(buf_color);

      SetUniformDomainSize(program);
      glUniform1ui(
          glGetUniformLocation(program, "pointSize"),
          kRadius * kInitHeight / 2);

      CHECK_ERROR();
      glDrawArrays(GL_POINTS, 0, particles.size());
      CHECK_ERROR();
      glUseProgram(0);
    };

    task2index["particles"] = tasks.size();
    tasks.push_back({program, render});
  }

  { // lines (frame, bonds, portals)
    GLuint program = CreateProgram(
        "../shaders/lines.vert.glsl", "../shaders/lines.frag.glsl",
        "../shaders/lines.geom.glsl");

    auto attr_point =
        std::make_shared<VertexAttribute<GLfloat, 2>>("point", program);
    auto attr_color =
        std::make_shared<VertexAttribute<GLfloat, 4>>("color", program);
    auto attr_width =
        std::make_shared<VertexAttribute<GLfloat, 1>>("width", program);

    auto render = [program, &ps = G->PS, attr_point, attr_color, attr_width]() {
      glUseProgram(program);

      std::vector<std::array<GLfloat, 2>> buf;
      std::vector<std::array<GLfloat, 4>> buf_color;
      std::vector<GLfloat> buf_width;
      size_t nprim = 0;

      using Color = std::array<Scal, 4>;
      auto add = [&](Vect a, Vect b, Color color, Scal width) {
        buf.push_back({a.x, a.y});
        buf_color.push_back(color);
        buf_width.push_back(width);

        buf.push_back({b.x, b.y});
        buf_color.push_back(color);
        buf_width.push_back(width);
        ++nprim;
      };
      auto rgb = [](float r, float g, float b) -> std::array<Scal, 4> {
        return {r, g, b, 1};
      };

      { // draw portals
        const auto& portals = ps->GetPortals();
        const Scal kPortalWidth = 0.015;
        const auto blue = rgb(0., .5, 1.);
        const auto orange = rgb(1., .5, 0.);
        for (auto& pair : portals) {
          add(pair[0].begin, pair[0].end, blue, kPortalWidth);
          add(pair[1].begin, pair[1].end, orange, kPortalWidth);
        }
        if (ps->portal_stage_ == 0) {
          if (ps->portal_mouse_moving_) {
            add(ps->portal_begin_, ps->portal_current_, blue, kPortalWidth);
          }
        } else {
          add(ps->portal_prev_.first, ps->portal_prev_.second, blue,
              kPortalWidth);
          if (ps->portal_mouse_moving_) {
            add(ps->portal_begin_, ps->portal_current_, orange, kPortalWidth);
          }
        }
      }

      { // draw bonds
        const Scal kBondWidth = 0.015;
        const auto& nr = ps->GetNoRendering();
        const auto& pos = ps->GetBlockData().position;
        const auto& bbi = ps->GetBlockById();
        for (auto bond : ps->GetBonds()) {
          const auto& a = bbi[bond.first];
          const auto& b = bbi[bond.second];
          if (!nr.count(bond)) {
            add(pos[a.first][a.second], pos[b.first][b.second], rgb(1, 1, 1),
                kBondWidth);
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
        const Scal kFrameWidth = 0.01;
        auto dom = ps->GetDomain();
        Vect dom0 = dom.A;
        Vect dom1 = dom.B;
        add(Vect(dom0.x, dom0.y), Vect(dom1.x, dom0.y), rgb(1, 1, 1),
            kFrameWidth);
        add(Vect(dom1.x, dom0.y), Vect(dom1.x, dom1.y), rgb(1, 1, 1),
            kFrameWidth);
        add(Vect(dom1.x, dom1.y), Vect(dom0.x, dom1.y), rgb(1, 1, 1),
            kFrameWidth);
        add(Vect(dom0.x, dom1.y), Vect(dom0.x, dom0.y), rgb(1, 1, 1),
            kFrameWidth);
      }

      fassert_equal(buf.size(), nprim * 2);
      fassert_equal(buf_width.size(), nprim * 2);
      fassert_equal(buf_color.size(), nprim * 2);

      attr_point->SetData(buf);
      attr_width->SetData(buf_width);
      attr_color->SetData(buf_color);

      SetUniformDomainSize(program);
      CHECK_ERROR();
      glDrawArrays(GL_LINES, 0, nprim * 2);
      CHECK_ERROR();
      glUseProgram(0);
    };

    task2index["lines"] = tasks.size();
    tasks.push_back({program, render});
  }

  { // blur rectangle
    GLuint program = CreateProgram(
        "../shaders/blur.vert.glsl", "../shaders/blur.frag.glsl", "");

    auto attr_point =
        std::make_shared<VertexAttribute<GLfloat, 2>>("point", program);

    GLuint tex = 0;
    glGenTextures(1, &tex);
    // const auto target = GL_TEXTURE_RECTANGLE_ARB;
    const auto target = GL_TEXTURE_2D;
    glBindTexture(target, tex);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_CLAMP);

    auto render = [program, attr_point, tex]() {
      glUseProgram(program);

      glBindTexture(target, tex);
      glReadBuffer(GL_FRONT);
      glCopyTexImage2D(target, 0, GL_RGB5, 0, 0, width, height, 0);

      const size_t nprim = 4;
      std::vector<std::array<GLfloat, 2>> buf{
          {-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
      attr_point->SetData(buf);

      glUniform2ui(glGetUniformLocation(program, "screenSize"), width, height);

      CHECK_ERROR();
      glDrawArrays(GL_TRIANGLE_FAN, 0, nprim);
      CHECK_ERROR();
      glUseProgram(0);
    };

    task2index["blur"] = tasks.size();
    tasks.push_back({program, render});
  }

  { // gui
    std::vector<std::uint32_t> colors_geo = {
        0xFF1F5B, 0x00CD6C, 0x009ADE, 0xAF58BA, 0xFFC61E, 0xF28522,
        0xA0B1BA, 0xA6761D, 0xE9002D, 0xFFAA00, 0x00B000};
    auto split = [](std::uint32_t c) -> std::array<float, 3> {
      std::array<float, 3> d;
      d[0] = ((c >> 16) & 0xFF) / 255.;
      d[1] = ((c >> 8) & 0xFF) / 255.;
      d[2] = ((c >> 0) & 0xFF) / 255.;
      return d; //
    };

    gui.AddButton(Gui::Button("abc", []() {}, {1, 1, 1}));
    gui.AddButton(Gui::Button("0abc", []() {}, split(colors_geo[0])));
    gui.AddButton(Gui::Button("1abc", []() {}, split(colors_geo[1])));
    gui.AddButton(Gui::Button("2abc", []() {}, split(colors_geo[2])));
    gui.SetWindowSize(width, height);

    auto font = LoadFont("../assets/font/font");
    auto tex_font = std::make_shared<Texture>();
    tex_font->SetData(font.data, font.width, font.height);

    GLuint program = CreateProgram(
        "../shaders/gui.vert.glsl", "../shaders/gui.frag.glsl",
        "../shaders/gui.geom.glsl");

    auto attr_lowcorner =
        std::make_shared<VertexAttribute<GLfloat, 2>>("lowcorner", program);
    auto attr_size =
        std::make_shared<VertexAttribute<GLfloat, 2>>("size", program);
    auto attr_color =
        std::make_shared<VertexAttribute<GLfloat, 4>>("color", program);
    auto attr_text =
        std::make_shared<VertexAttribute<GLchar, 4>>("text", program);

    auto render = [program, &ps = G->PS, attr_lowcorner, attr_color, attr_size,
                   attr_text, tex_font, font]() {
      glUseProgram(program);

      std::vector<std::array<GLfloat, 2>> buf_lowcorner;
      std::vector<std::array<GLfloat, 2>> buf_size;
      std::vector<std::array<GLfloat, 4>> buf_color;
      std::vector<std::array<GLchar, 4>> buf_text;
      size_t nprim = 0;

      using Color = std::array<Scal, 4>;
      auto add = [&](Vect lowcorner, Vect size, Color color, std::string text) {
        buf_lowcorner.push_back(lowcorner);
        buf_size.push_back(size);
        buf_color.push_back(color);
        typename decltype(buf_text)::value_type q;
        for (size_t i = 0; i < q.size(); ++i) {
          if (i < text.size()) {
            q[i] = text[i] - 0x20;
          } else {
            q[i] = ' ' - 0x20;
          }
        }
        buf_text.push_back(q);
        ++nprim;
      };
      auto rgba = [](std::array<float, 3> c) -> std::array<Scal, 4> {
        return {c[0], c[1], c[2], 1};
      };

      auto& buttons = gui.GetButtons();

      for (auto& b : buttons) {
        add(b.lowcorner, b.size, rgba(b.color), b.label);
      }

      fassert_equal(buf_lowcorner.size(), nprim);
      fassert_equal(buf_size.size(), nprim);
      fassert_equal(buf_color.size(), nprim);
      fassert_equal(buf_text.size(), nprim);

      tex_font->Bind();

      attr_lowcorner->SetData(buf_lowcorner);
      attr_size->SetData(buf_size);
      attr_color->SetData(buf_color);
      attr_text->SetDataInt(buf_text);

      glUniform2ui(
          glGetUniformLocation(program, "font_size"), font.width, font.height);
      glUniform1f(glGetUniformLocation(program, "font_offset"), font.offset);

      SetUniformDomainSize(program);
      CHECK_ERROR();
      glDrawArrays(GL_POINTS, 0, nprim);
      CHECK_ERROR();
      glUseProgram(0);
    };

    task2index["gui"] = tasks.size();
    tasks.push_back({program, render});
  }

  tasks_indices = {
      task2index["blur"],
      task2index["particles"],
      task2index["lines"],
      task2index["gui"],
  };

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
            gui.SetWindowSize(width, height);
            glViewport(0, 0, width, height);
            G->SetWindowSize(width, height);
            glDrawBuffer(GL_FRONT);
            glClear(GL_COLOR_BUFFER_BIT);
            glDrawBuffer(GL_BACK);
            break;
        }
      }
    }

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
