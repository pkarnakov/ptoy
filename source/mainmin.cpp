#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <GL/glew.h>
#include <SDL.h>
#include <SDL_opengl.h>

#include "game.h"
#include "logger.h"
#include "view_gl.h"

enum class MouseState {
  None,
  Attraction,
  Repulsion,
  Bond,
  Pick,
  Freeze,
  Portal,
};

const char* MouseStateName(MouseState s) {
  switch (s) {
    case MouseState::None:
      return "none";
    case MouseState::Attraction:
      return "attraction";
    case MouseState::Repulsion:
      return "repulsion";
    case MouseState::Bond:
      return "bonds";
    case MouseState::Pick:
      return "pick";
    case MouseState::Freeze:
      return "freeze";
    case MouseState::Portal:
      return "portal";
    default:
      return nullptr;
  }
}

class Gui {
 public:
  using Color = std::array<float, 3>;

  static constexpr unsigned kButtonHeight = 40;
  static constexpr unsigned kButtonWidth = 40;

  struct Window {
    Vect size;
  };

  struct Button {
    std::string label;
    std::function<void()> handler;
    Color color;
    Vect lowcorner;
    Vect size;

    Button(std::string label_, std::function<void()> handler_, Color color_)
        : label(label_), handler(handler_), color(color_) {
      size = Vect(kButtonWidth, kButtonHeight);
    }
  };

  Button* FindButton(int x, int y) {
    for (auto& b : buttons_) {
      if (b.lowcorner.x < x && x <= b.lowcorner.x + b.size.x &&
          b.lowcorner.y < y && y <= b.lowcorner.y + b.size.y) {
        return &b;
      }
    }
    return nullptr;
  }

  void AddButton(const Button& e) {
    buttons_.push_back(e);
  }
  std::list<Button>& GetButtons() {
    return buttons_;
  }
  const std::list<Button>& GetButtons() const {
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

  std::list<Button> buttons_;
  Window window_;
};

std::chrono::milliseconds last_frame_time;
Scal last_frame_game_time;
std::chrono::milliseconds last_report_time;
std::atomic<Scal> next_game_time_target;

std::unique_ptr<Game> gameinst;

Gui gui;

std::map<MouseState, Gui::Button*> mouse_to_button;

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
  std::chrono::milliseconds current_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now().time_since_epoch());
  std::chrono::milliseconds time_past_from_last_frame =
      current_time - last_frame_time;

  auto new_frame_time = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now().time_since_epoch());
  auto new_frame_game_time = gameinst->partsys->GetTime();

  const Scal frame_real_duration_s =
      (new_frame_time - last_frame_time).count() / 1000.;

  if ((new_frame_time - last_report_time).count() > 1000.) {
    std::cout << "fps: " << 1. / frame_real_duration_s << ", game rate="
              << (new_frame_game_time - last_frame_game_time) /
                     frame_real_duration_s
              << ", particles="
              << gameinst->partsys->GetNumParticles()
              //<< ", max_per_cell="
              //<< gameinst->partsys->GetNumPerCell()
              << ", t=" << gameinst->partsys->GetTime()
              << ", steps=" << gameinst->partsys->GetNumSteps() << std::endl;
    last_report_time = new_frame_time;
    // gameinst->partsys->Blocks.print_status();
  }

  const Scal game_rate_target = 1.5;
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
    std::lock_guard<std::mutex> lg(gameinst->partsys->m_buffer_);
    gameinst->partsys->SetRendererReadyForNext(true);
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
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_nested(0);
  // omp_set_num_threads(std::thread::hardware_concurrency());
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
    gameinst->partsys->step(next_game_time_target, pause);
    if (pause) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
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

  gameinst = std::unique_ptr<Game>(new Game(width, height));

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

  glewInit();

  std::vector<std::uint32_t> colors_geo = {
      0xFF1F5B, 0x00CD6C, 0x009ADE, 0xAF58BA, 0xFFC61E, 0xF28522,
      0xA0B1BA, 0xA6761D, 0xE9002D, 0xFFAA00, 0x00B000};

  auto split_rgb = [](std::uint32_t c) -> std::array<float, 3> {
    std::array<float, 3> d;
    d[0] = ((c >> 16) & 0xFF) / 255.;
    d[1] = ((c >> 8) & 0xFF) / 255.;
    d[2] = ((c >> 0) & 0xFF) / 255.;
    return d;
  };

  auto gray = split_rgb(colors_geo[6]);
  glClearColor(gray[0] * 0.5, gray[1] * 0.5, gray[2] * 0.5, 1.0);
  glDrawBuffer(GL_FRONT);
  glClear(GL_COLOR_BUFFER_BIT);
  glDrawBuffer(GL_BACK);

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glEnable(GL_PROGRAM_POINT_SIZE);

  glProvokingVertex(GL_FIRST_VERTEX_CONVENTION);

  std::map<std::string, size_t> task2index;

  MouseState mouse_state = MouseState::Repulsion;
  bool help_state = true;

  auto set_mouse_mode = [&](MouseState s) {
    mouse_state = s;
    std::cout << std::string("Mouse switched to mode: ") + MouseStateName(s)
              << std::endl;
    switch (s) {
      case MouseState::Repulsion:
        gameinst->partsys->SetForceAttractive(false);
        break;
      case MouseState::Attraction:
        gameinst->partsys->SetForceAttractive(true);
        break;
      default:
        break;
    }
  };

  { // Draw particles.
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

    auto render = [&, program, tex_circle, &ps = gameinst->partsys, attr_point,
                   attr_color]() {
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
      auto rgb = split_rgb(colors_geo[1]);
      glUniform3f(
          glGetUniformLocation(program, "facecolor"), rgb[0], rgb[1], rgb[2]);

      CHECK_ERROR();
      glDrawArrays(GL_POINTS, 0, particles.size());
      CHECK_ERROR();
      glUseProgram(0);
    };

    task2index["particles"] = tasks.size();
    tasks.push_back({program, render});
  }

  { // Draw lines (frame, bonds, portals).
    GLuint program = CreateProgram(
        "../shaders/lines.vert.glsl", "../shaders/lines.frag.glsl",
        "../shaders/lines.geom.glsl");

    auto attr_point =
        std::make_shared<VertexAttribute<GLfloat, 2>>("point", program);
    auto attr_color =
        std::make_shared<VertexAttribute<GLfloat, 4>>("color", program);
    auto attr_width =
        std::make_shared<VertexAttribute<GLfloat, 1>>("width", program);

    auto render = [&, program, &ps = gameinst->partsys, attr_point, attr_color,
                   attr_width]() {
      glUseProgram(program);

      std::vector<std::array<GLfloat, 2>> buf;
      std::vector<std::array<GLfloat, 4>> buf_color;
      std::vector<GLfloat> buf_width;
      size_t nprim = 0;

      using Color = std::array<Scal, 4>;
      auto add = [&](Vect a, Vect b, Color c, Scal w) {
        buf.push_back({a.x, a.y});
        buf_color.push_back(c);
        buf_width.push_back(w);

        buf.push_back({b.x, b.y});
        buf_color.push_back(c);
        buf_width.push_back(w);
        ++nprim;
      };
      auto rgb = [](float r, float g, float b) -> Color {
        return {r, g, b, 1};
      };
      auto rgba = [](std::array<Scal, 3> c) -> Color {
        return {c[0], c[1], c[2], 1};
      };

      { // Draw portals.
        const auto& portals = ps->GetPortals();
        const Scal kPortalWidth = 0.015;
        const auto blue = rgba(split_rgb(colors_geo[2]));
        const auto orange = rgba(split_rgb(colors_geo[5]));
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

      { // Draw bonds.
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

      { // Draw frozen particles.
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

      { // Draw frame.
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

  { // Blur rectangle.
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

  { // Initialize GUI.
    auto add_button = [&](std::string lbl, MouseState s) {
      gui.AddButton(Gui::Button( //
          lbl, [&, s]() { set_mouse_mode(s); }, {1, 1, 1}));
      mouse_to_button[s] = &gui.GetButtons().back();
    };
    add_button("R", MouseState::Repulsion);
    add_button("A", MouseState::Attraction);
    add_button("P", MouseState::Pick);
    add_button("F", MouseState::Freeze);
    add_button("O", MouseState::Portal);
    add_button("B", MouseState::Bond);

    gui.SetWindowSize(width, height);
  }

  { // GUI.
    GLuint program = CreateProgram(
        "../shaders/gui.vert.glsl", "../shaders/gui.frag.glsl",
        "../shaders/gui.geom.glsl");

    auto attr_lowcorner =
        std::make_shared<VertexAttribute<GLfloat, 2>>("lowcorner", program);
    auto attr_size =
        std::make_shared<VertexAttribute<GLfloat, 2>>("size", program);
    auto attr_color =
        std::make_shared<VertexAttribute<GLfloat, 4>>("color", program);

    auto render = [&, program, &ps = gameinst->partsys, attr_lowcorner,
                   attr_color, attr_size]() {
      glUseProgram(program);

      std::vector<std::array<GLfloat, 2>> buf_lowcorner;
      std::vector<std::array<GLfloat, 2>> buf_scale;
      std::vector<std::array<GLfloat, 4>> buf_color;
      size_t nprim = 0;

      using Color = std::array<Scal, 4>;
      auto add = [&](Vect lowcorner, Vect size, Color color) {
        buf_lowcorner.push_back(lowcorner);
        buf_scale.push_back(size);
        buf_color.push_back(color);
        ++nprim;
      };
      auto rgba = [](std::array<float, 3> c) -> std::array<Scal, 4> {
        return {c[0], c[1], c[2], 1};
      };

      auto& buttons = gui.GetButtons();

      for (auto& b : buttons) {
        b.color = {1., 1., 1.};
      }
      if (mouse_to_button.count(mouse_state)) {
        mouse_to_button[mouse_state]->color = split_rgb(colors_geo[0]);
      }

      for (auto& b : buttons) {
        add(b.lowcorner, b.size, rgba(b.color));
      }

      attr_lowcorner->SetData(buf_lowcorner);
      attr_size->SetData(buf_scale);
      attr_color->SetData(buf_color);

      SetUniformDomainSize(program);
      CHECK_ERROR();
      glDrawArrays(GL_POINTS, 0, nprim);
      CHECK_ERROR();
      glUseProgram(0);
    };

    task2index["gui"] = tasks.size();
    tasks.push_back({program, render});
  }

  { // Text.
    auto font = LoadFont("../assets/font/font");
    auto tex_font = std::make_shared<Texture>();
    tex_font->SetData(font.mask, font.width, font.height);

    GLuint program = CreateProgram(
        "../shaders/text.vert.glsl", "../shaders/text.frag.glsl",
        "../shaders/text.geom.glsl");

    auto attr_lowcorner =
        std::make_shared<VertexAttribute<GLfloat, 2>>("lowcorner", program);
    auto attr_size =
        std::make_shared<VertexAttribute<GLfloat, 1>>("size", program);
    auto attr_color =
        std::make_shared<VertexAttribute<GLfloat, 4>>("color", program);
    auto attr_char =
        std::make_shared<VertexAttribute<GLfloat, 1>>("char", program);
    auto attr_width =
        std::make_shared<VertexAttribute<GLfloat, 1>>("width", program);

    auto render = [program, &ps = gameinst->partsys, attr_lowcorner, attr_color,
                   attr_size, attr_char, attr_width, tex_font, font,
                   &mouse_state, &help_state]() {
      glUseProgram(program);

      std::vector<std::array<GLfloat, 2>> buf_lowcorner;
      std::vector<GLfloat> buf_scale;
      std::vector<std::array<GLfloat, 4>> buf_color;
      std::vector<GLfloat> buf_char;
      std::vector<GLfloat> buf_width;
      size_t nprim = 0;

      using Color = std::array<Scal, 4>;
      auto add_text = [&](Vect lowcorner, Scal scale, Color color,
                          std::string text) {
        Scal x = 0;
        Scal y = 0;
        for (auto c : text) {
          if (c == '\n') {
            x = 0;
            y -= font.height * 1.25 * scale;
          } else {
            buf_lowcorner.push_back(lowcorner + Vect(x, y));
            buf_scale.push_back(scale);
            buf_color.push_back(color);
            buf_char.push_back(font.xoffsets[c]);
            buf_width.push_back(font.widths[c]);
            ++nprim;
            x += font.widths[c] * scale;
          }
        }
      };
      auto rgba = [](std::array<float, 3> c) -> std::array<Scal, 4> {
        return {c[0], c[1], c[2], 1};
      };

      auto& buttons = gui.GetButtons();
      for (auto& b : buttons) {
        add_text(
            b.lowcorner + (b.size - font.GetLength(b.label)) * 0.5, 1,
            rgba({0, 0, 0}), b.label);
      }

      if (help_state) {
        add_text(
            buttons.front().lowcorner + Vect(5, -50), 0.5, rgba({1, 1, 1}),
            R"(R: switch to repulsion
A: switch to attraction
P: switch to pick
F: switch to freeze
O: switch to portal
B: switch to bond
N: switch to no action
G: toggle gravity
I: remove last pair of portals
SPACE: toggle pause
H: show this help
Q: quit after three presses
          )");
      }

      tex_font->Bind();

      attr_lowcorner->SetData(buf_lowcorner);
      attr_size->SetData(buf_scale);
      attr_color->SetData(buf_color);
      attr_char->SetData(buf_char);
      attr_width->SetData(buf_width);

      glUniform2ui(
          glGetUniformLocation(program, "font_size"), font.width, font.height);

      SetUniformDomainSize(program);
      CHECK_ERROR();
      glDrawArrays(GL_POINTS, 0, nprim);
      CHECK_ERROR();
      glUseProgram(0);
    };

    task2index["text"] = tasks.size();
    tasks.push_back({program, render});
  }

  tasks_indices = {
      // task2index["blur"],
      task2index["particles"],
      task2index["lines"],
      task2index["gui"],
      task2index["text"],
  };

  std::thread computation_thread(cycle);

  gameinst->partsys->SetForce(Vect(0., 0.), false);

  // Main loop flag.
  quit = false;

  // Event handler.
  SDL_Event e;

  size_t quit_count = 0;

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
          case 'h':
            help_state = !help_state;
            std::cout << "Toggle help: " << (help_state ? "on" : "off")
                      << std::endl;
            break;
          case 'n':
            set_mouse_mode(MouseState::None);
            break;
          case 'r':
            set_mouse_mode(MouseState::Repulsion);
            break;
          case 'f':
            set_mouse_mode(MouseState::Freeze);
            break;
          case 'p':
            set_mouse_mode(MouseState::Pick);
            break;
          case 'a':
            set_mouse_mode(MouseState::Attraction);
            break;
          case 'b':
            set_mouse_mode(MouseState::Bond);
            break;
          case 'o':
            set_mouse_mode(MouseState::Portal);
            break;
          case 'i':
            std::cout << "Remove last portal" << std::endl;
            gameinst->partsys->RemoveLastPortal();
            break;
          case 'g':
            gameinst->partsys->InvertGravity();
            std::cout << (gameinst->partsys->GetGravity() ? "Gravity on"
                                                          : "Gravity off")
                      << std::endl;
            break;
          case ' ':
            pause = !pause;
            break;
        }
      } else if (e.type == SDL_MOUSEMOTION) {
        switch (mouse_state) {
          case MouseState::Attraction:
          case MouseState::Repulsion:
            gameinst->partsys->SetForce(GetDomainMousePosition());
            break;
          case MouseState::Bond:
            gameinst->partsys->BondsMove(GetDomainMousePosition());
            break;
          case MouseState::Pick:
            gameinst->partsys->PickMove(GetDomainMousePosition());
            break;
          case MouseState::Freeze:
            gameinst->partsys->FreezeMove(GetDomainMousePosition());
            break;
          case MouseState::Portal:
            gameinst->partsys->PortalMove(GetDomainMousePosition());
            break;
          case MouseState::None:
            break;
        }
      } else if (e.type == SDL_MOUSEBUTTONDOWN) {
        int x, y;
        SDL_GetMouseState(&x, &y);
        if (auto b = gui.FindButton(x, height - y)) {
          std::cout << "found button at " << NAMEVALUE(x) << " " << NAMEVALUE(y)
                    << " with " << NAMEVALUE(b->label) << std::endl;
          if (b->handler) {
            b->handler();
          }
        } else {
          switch (mouse_state) {
            case MouseState::Attraction:
            case MouseState::Repulsion:
              gameinst->partsys->SetForce(GetDomainMousePosition(), true);
              break;
            case MouseState::Bond:
              gameinst->partsys->BondsStart(GetDomainMousePosition());
              break;
            case MouseState::Pick:
              gameinst->partsys->PickStart(GetDomainMousePosition());
              break;
            case MouseState::Freeze:
              gameinst->partsys->FreezeStart(GetDomainMousePosition());
              break;
            case MouseState::Portal:
              gameinst->partsys->PortalStart(GetDomainMousePosition());
              break;
            case MouseState::None:
              break;
          }
        }
      } else if (e.type == SDL_MOUSEBUTTONUP) {
        switch (mouse_state) {
          case MouseState::Attraction:
          case MouseState::Repulsion:
            gameinst->partsys->SetForce(false);
            break;
          case MouseState::Bond:
            gameinst->partsys->BondsStop(GetDomainMousePosition());
            break;
          case MouseState::Pick:
            gameinst->partsys->PickStop(GetDomainMousePosition());
            break;
          case MouseState::Freeze:
            gameinst->partsys->FreezeStop(GetDomainMousePosition());
            break;
          case MouseState::Portal:
            gameinst->partsys->PortalStop(GetDomainMousePosition());
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
            gameinst->SetWindowSize(width, height);
            glDrawBuffer(GL_FRONT);
            glClear(GL_COLOR_BUFFER_BIT);
            glDrawBuffer(GL_BACK);
            break;
        }
      }
    }

    glClear(GL_COLOR_BUFFER_BIT);
    display();
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    SDL_GL_SwapWindow(window);
  }

  // Finalize
  SDL_GL_DeleteContext(glcontext);
  SDL_DestroyWindow(window);
  SDL_Quit();

  computation_thread.join();

  return 0;
}
