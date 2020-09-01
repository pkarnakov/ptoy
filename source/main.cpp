#include <iostream>

#include <GL/glew.h>
#include <SDL.h>
#include <SDL_opengl.h>
#include "game.hpp"
#include "logger.h"

#include <atomic>
#include <chrono>
#include <thread>
#include <fstream>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

// TODO: reorder header includes (std first)
// TODO: fix formatting

using namespace std::chrono;

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
} // namespace wrap

milliseconds last_frame_time;
Scal last_frame_game_time;
milliseconds last_report_time;
std::atomic<Scal> next_game_time_target;

std::unique_ptr<game> G;

GLdouble width, height; /* window width and height */

std::atomic<bool> flag_display;
std::atomic<bool> quit;
std::atomic<bool> pause;

using std::cout;
using std::endl;

GLuint gProgramID = 0;
GLint gPointArray = -1;
GLint gColorArray = -1;
GLint gScreenSizeLocation = -1;
GLuint gVBO_point = 0;
GLuint gVBO_color = 0;

int frame_number;

void init() {
  width = 800.0; /* initial window width and height, */
  height = 800.0; /* within which we draw. */

  last_frame_game_time = 0.;
  pause = false;
}

/* Callback functions for GLUT */

/* Draw the window - this is where all the GL actions are */
void display(void) {
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

  /* clear the screen to white */
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(gProgramID);

  auto& particles = G->PS->GetParticles();
  std::vector<GLfloat> buf(particles.size() * 2);
  std::vector<GLfloat> buf_color(particles.size());

  {
    std::lock_guard<std::mutex> lg(G->PS->m_buffer_);
    for (size_t i = 0; i < particles.size(); ++i) {
      auto& p = particles[i];
      buf[i * 2] = p.p.x;
      buf[i * 2 + 1] = p.p.y;
      Scal f = 0.5 + p.v.length() / 7.; // color intensity
      f = std::min<Scal>(std::max<Scal>(f, 0.), 1.);
      buf_color[i] = f;
    }
  }

  glBindBuffer(GL_ARRAY_BUFFER, gVBO_point);
  GLint size = 0;
  glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
  if (size < int(buf.size())) { // reallocate
    glBindBuffer(GL_ARRAY_BUFFER, gVBO_point);
    glBufferData(
        GL_ARRAY_BUFFER, buf.size() * sizeof(GLfloat), buf.data(),
        GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, gVBO_color);
    glBufferData(
        GL_ARRAY_BUFFER, buf_color.size() * sizeof(GLfloat), buf_color.data(),
        GL_DYNAMIC_DRAW);
  } else { // use existing
    glBindBuffer(GL_ARRAY_BUFFER, gVBO_point);
    glBufferSubData(
        GL_ARRAY_BUFFER, 0, buf.size() * sizeof(GLfloat), buf.data());
    glBindBuffer(GL_ARRAY_BUFFER, gVBO_color);
    glBufferSubData(
        GL_ARRAY_BUFFER, 0, buf_color.size() * sizeof(GLfloat),
        buf_color.data());
  }

  glEnableVertexAttribArray(gPointArray);
  glEnableVertexAttribArray(gColorArray);

  glBindBuffer(GL_ARRAY_BUFFER, gVBO_point);
  glVertexAttribPointer(
      gPointArray, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), NULL);
  glBindBuffer(GL_ARRAY_BUFFER, gVBO_color);
  glVertexAttribPointer(
      gColorArray, 1, GL_FLOAT, GL_FALSE, 1 * sizeof(GLfloat), NULL);

  vect A(-1., -1.), B(-1 + 2. * width / 800, -1. + 2. * height / 800);
  glUniform2ui(glGetUniformLocation(gProgramID, "screenSize"), width, height);
  glUniform1ui(
      glGetUniformLocation(gProgramID, "pointSize"), kRadius * 800 * 0.5);
  glUniform2f(glGetUniformLocation(gProgramID, "domain0"), A.x, A.y);
  glUniform2f(glGetUniformLocation(gProgramID, "domain1"), B.x, B.y);

  glDrawArrays(GL_POINTS, 0, particles.size());

  glDisableVertexAttribArray(gPointArray);
  glDisableVertexAttribArray(gColorArray);

  glUseProgram(0);

  {
    std::lock_guard<std::mutex> lg(G->PS->m_buffer_);
    G->R->DrawAll();
    G->PS->SetRendererReadyForNext(true);
  }

  glFlush();

  flag_display = false;
}

void cycle() {
  // omp_set_dynamic(0);
  // omp_set_nested(0);
  // omp_set_num_threads(std::thread::hardware_concurrency());
#ifdef _OPENMP
  omp_set_num_threads(1);
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

vect GetDomainMousePosition(int x, int y) {
  vect c;
  vect A(-1., -1.), B(-1 + 2. * width / 800, -1. + 2. * height / 800);
  c.x = A.x + (B.x - A.x) * (x / width);
  c.y = B.y + (A.y - B.y) * (y / height);
  return c;
}

vect GetDomainMousePosition() {
  int x, y;
  SDL_GetMouseState(&x, &y);
  return GetDomainMousePosition(x, y);
}

void printShaderLog(GLuint shader) {
  // Make sure name is shader
  if (glIsShader(shader)) {
    // Shader log length
    int infoLogLength = 0;
    int maxLength = infoLogLength;

    // Get info string length
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);

    // Allocate string
    char* infoLog = new char[maxLength];

    // Get info log
    glGetShaderInfoLog(shader, maxLength, &infoLogLength, infoLog);
    if (infoLogLength > 0) {
      // Print Log
      printf("%s\n", infoLog);
    }

    // Deallocate string
    delete[] infoLog;
  } else {
    printf("Name %d is not a shader\n", shader);
  }
}
void printProgramLog(GLuint program) {
  // Make sure name is shader
  if (glIsProgram(program)) {
    // Program log length
    int infoLogLength = 0;
    int maxLength = infoLogLength;

    // Get info string length
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);

    // Allocate string
    char* infoLog = new char[maxLength];

    // Get info log
    glGetProgramInfoLog(program, maxLength, &infoLogLength, infoLog);
    if (infoLogLength > 0) {
      // Print Log
      printf("%s\n", infoLog);
    }

    // Deallocate string
    delete[] infoLog;
  } else {
    printf("Name %d is not a program\n", program);
  }
}

int main() {
  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    SDL_Log("Unable to initialize SDL: %s\n", SDL_GetError());
    return 1;
  }

  init();

  SDL_Window* window = SDL_CreateWindow(
      "ptoy", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height,
      SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

  if (window == NULL) {
    SDL_Log("Unable to create window: %s\n", SDL_GetError());
    return 1;
  }

  SDL_GLContext glcontext = SDL_GL_CreateContext(window);

  auto gray = 0.5;
  glClearColor(gray, gray, gray, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  SDL_GL_SwapWindow(window);

  glewInit();

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glEnable(GL_PROGRAM_POINT_SIZE);

  auto readfile = [](std::string path) -> std::string {
    std::ifstream fs(path);
    fassert(fs.good(), "can't open file + '" + path + "'");
    std::stringstream ss;
    ss << fs.rdbuf();
    return ss.str();
  };

  auto shd = [&] {
    gProgramID = glCreateProgram();
    // Create vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);

    const std::string src_vert = readfile("../shaders/s.vert");
    const GLchar* vertexShaderSource[] = {src_vert.c_str()};

    // Set vertex source
    glShaderSource(vertexShader, 1, vertexShaderSource, NULL);

    // Compile vertex source
    glCompileShader(vertexShader);

    // Check vertex shader for errors
    GLint vShaderCompiled = GL_FALSE;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &vShaderCompiled);
    if (vShaderCompiled != GL_TRUE) {
      printf("Unable to compile vertex shader %d!\n", vertexShader);
      printShaderLog(vertexShader);
      fassert(false);
    }

    // Attach vertex shader to program
    glAttachShader(gProgramID, vertexShader);

    // Create fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    // Get fragment source
    const std::string src_frag = readfile("../shaders/s.frag");
    const GLchar* fragmentShaderSource[] = {src_frag.c_str()};

    // Set fragment source
    glShaderSource(fragmentShader, 1, fragmentShaderSource, NULL);

    // Compile fragment source
    glCompileShader(fragmentShader);

    // Check fragment shader for errors
    GLint fShaderCompiled = GL_FALSE;
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &fShaderCompiled);
    if (fShaderCompiled != GL_TRUE) {
      printf("Unable to compile fragment shader %d!\n", fragmentShader);
      printShaderLog(fragmentShader);
      fassert(false);
    }
    glAttachShader(gProgramID, fragmentShader);
    glLinkProgram(gProgramID);
    GLint programSuccess = GL_TRUE;
    glGetProgramiv(gProgramID, GL_LINK_STATUS, &programSuccess);
    if (programSuccess != GL_TRUE) {
      printProgramLog(gProgramID);
      printf("Error linking program %d!\n", gProgramID);
      fassert(false);
    }
    gPointArray = wrap::glGetAttribLocation(gProgramID, "point");
    gColorArray = wrap::glGetAttribLocation(gProgramID, "color");
    glGenBuffers(1, &gVBO_point);
    glGenBuffers(1, &gVBO_color);
  };
  shd();

  G = std::unique_ptr<game>(new game(width, height));
  std::thread computation_thread(cycle);

  G->PS->SetForce(vect(0., 0.), false);

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

    display();
    SDL_GL_SwapWindow(window);
    SDL_Delay(30);
  }

  // Finalize
  SDL_GL_DeleteContext(glcontext);
  SDL_DestroyWindow(window);
  SDL_Quit();

  computation_thread.join();

  return 0;
}
