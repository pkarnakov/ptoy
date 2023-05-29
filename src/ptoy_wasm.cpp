#include <stdint.h>
#include <stdlib.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>

#include "control.h"
#include "game.h"
#include "logger.h"
#include "macros.h"
#include "scene.h"

struct SceneData {
  struct Particles {
    std::vector<Vect> p;
    std::vector<Vect> v;
  };
  Particles particles;

  using Portal = Scene::Portal;
  std::vector<std::array<Portal, 2>> portals;
};

const int g_width = 800;
const int g_height = 800;
Scene g_scene;
SceneData g_data;
std::shared_ptr<Game> g_gameinst;
std::shared_ptr<Control> g_control;
std::string g_buf;
bool state_pause;

void UpdateScene() {
  auto gameinst = g_gameinst;
  const auto& particles = gameinst->partsys->GetParticles();
  g_data.particles.p.resize(particles.size());
  g_data.particles.v.resize(particles.size());
  for (size_t i = 0; i < particles.size(); ++i) {
    g_data.particles.p[i] = particles[i].p;
    g_data.particles.v[i] = particles[i].v;
  }
  g_scene.particles.p = g_data.particles.p;
  g_scene.particles.v = g_data.particles.v;

  const auto& portals = gameinst->partsys->GetPortals();
  g_data.portals.resize(portals.size());
  for (size_t i = 0; i < portals.size(); ++i) {
    for (size_t j : {0, 1}) {
      g_data.portals[i][j].pa = portals[i][j].begin;
      g_data.portals[i][j].pb = portals[i][j].end;
    }
  }
  g_scene.portals = g_data.portals;
}

static void main_loop() {
  auto gameinst = g_gameinst;
  gameinst->partsys->SetRendererReadyForNext(true);
  const auto dt = 0.01;
  gameinst->partsys->step(gameinst->partsys->GetTime() + dt, state_pause);
  UpdateScene();
  EM_ASM_({ draw(); });
}

extern "C" {
int SetConfig(const char*) {
  return 0;
}
const char* GetConfig() {
  return "";
}
int GetParticles(uint16_t* data, int max_size) {
  auto gameinst = g_gameinst;
  const int entrysize = 2;
  int i = 0;
  for (auto p : g_scene.particles.p) {
    if (i + entrysize > max_size) {
      break;
    }
    data[i + 0] = (1 + p.x) * gameinst->width_ * 0.5;
    data[i + 1] = (1 - p.y) * gameinst->height_ * 0.5;
    i += entrysize;
  }
  return i;
}
int GetPortals(uint16_t* data, int max_size) {
  auto gameinst = g_gameinst;
  const int entrysize = 8;
  int i = 0;
  for (auto pair : g_scene.portals) {
    if (i + entrysize > max_size) {
      break;
    }
    auto append = [&](Vect p) {
      data[i + 0] = (1 + p.x) * gameinst->width_ * 0.5;
      data[i + 1] = (1 - p.y) * gameinst->height_ * 0.5;
      i += 2;
    };
    append(pair[0].pa);
    append(pair[0].pb);
    append(pair[1].pa);
    append(pair[1].pb);
  }
  return i;
}
void SendKeyDown(char keysym) {
  auto control = g_control;
  control->SendKeyDown(keysym);
}
void SendMouseMotion(float x, float y) {
  auto control = g_control;
  control->SendMouseMotion({x, y});
}
void SendMouseDown(float x, float y) {
  auto control = g_control;
  control->SendMouseDown({x, y});
}
void SendMouseUp(float x, float y) {
  auto control = g_control;
  control->SendMouseUp({x, y});
}
void SetControlDebug(int flag) {
  auto control = g_control;
  control->debug = flag;
}
void Init() {
  state_pause = false;
  g_gameinst = std::make_shared<Game>(g_width, g_height);
  g_control = std::make_shared<Control>(g_gameinst->partsys.get());
  UpdateScene();
}
void SetPause(int flag) {
  state_pause = flag;
}
int GetGravity() {
  auto gameinst = g_gameinst;
  return gameinst->partsys->GetGravity();
}
void SetGravity(int flag) {
  auto gameinst = g_gameinst;
  gameinst->partsys->SetGravity(flag);
}
const char* GetMouseMode() {
  auto gameinst = g_gameinst;
  auto control = g_control;
  return Control::MouseModeToStr(control->mouse_mode);
}
} // extern "C"

int main() {
  Init();
  emscripten_set_canvas_element_size("#canvas", g_width, g_height);
  emscripten_set_main_loop(main_loop, 30, 1);
}
