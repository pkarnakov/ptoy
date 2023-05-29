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

static constexpr int kScale = 1;

struct SceneData {
  struct Particles {
    std::vector<Vect> p;
    std::vector<Vect> v;
  };
  Particles particles;

  using Portal = Scene::Portal;
  std::vector<std::array<Portal, 2>> portals;
};

Scene g_scene;
SceneData g_data;
std::unique_ptr<Game> gameinst;
std::unique_ptr<Control> control;
std::string g_buf;
bool state_pause;
bool state_quit;

void UpdateScene() {
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
  gameinst->partsys->SetRendererReadyForNext(true);
  const auto dt = 0.02;
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
  const int entrysize = 2;
  int i = 0;
  for (auto p : g_scene.particles.p) {
    if (i + entrysize > max_size) {
      break;
    }
    data[i + 0] = (1 + p.x) * gameinst->width_ * 0.5 / kScale;
    data[i + 1] = (1 - p.y) * gameinst->height_ * 0.5 / kScale;
    i += entrysize;
  }
  return i;
}
int GetPortals(uint16_t* data, int max_size) {
  const int entrysize = 8;
  int i = 0;
  for (auto pair : g_scene.portals) {
    if (i + entrysize > max_size) {
      break;
    }
    auto append = [&](Vect p) {
      data[i + 0] = (1 + p.x) * gameinst->width_ * 0.5 / kScale;
      data[i + 1] = (1 - p.y) * gameinst->height_ * 0.5 / kScale;
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
  control->SendKeyDown(keysym);
}
void SendMouseMotion(float x, float y) {
  control->SendMouseMotion({x, y});
}
void SendMouseDown(float x, float y) {
  control->SendMouseDown({x, y});
}
void SendMouseUp(float x, float y) {
  control->SendMouseUp({x, y});
}
} // extern "C"

int main() {
  state_pause = false;
  state_quit = false;
  const int width = 800;
  const int height = 800;
  gameinst = std::make_unique<Game>(width, height);
  control = std::make_unique<Control>(gameinst->partsys.get());
  //control->debug = true;
  UpdateScene();

  emscripten_set_canvas_element_size(
      "#canvas", width / kScale, height / kScale);
  emscripten_set_main_loop(main_loop, 30, 1);
}
