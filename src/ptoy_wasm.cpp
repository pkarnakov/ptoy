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

#include "game.h"
#include "logger.h"
#include "macros.h"
#include "scene.h"
#include "view_text.h"

static constexpr int kScale = 2;

struct SceneData {
  struct Particles {
    std::vector<Vect> p;
    std::vector<Vect> v;
  };
  Particles particles;

  using Portal = Scene::Portal;
  std::vector<std::array<Portal, 2>> portals;
};

std::unique_ptr<View> g_view;
Scene g_scene;
SceneData g_data;
std::unique_ptr<Game> gameinst;
std::string g_buf;
std::atomic<bool> state_pause;
std::atomic<bool> state_quit;

static void main_loop() {
  gameinst->partsys->SetRendererReadyForNext(true);
  const auto dt = 0.02;
  gameinst->partsys->step(gameinst->partsys->GetTime() + dt, state_pause);
  EM_ASM_({ Draw(); });
}

extern "C" {
int SetConfig(const char*) {
  return 0;
}

const char* GetConfig() {
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
  g_view->SetScene(g_scene);
  std::stringstream buf;
  buf << gameinst->partsys->GetTime() << '\n';
  for (auto p : g_scene.particles.p) {
    buf << p;
    break;
  }
  g_buf = buf.str();
  return g_buf.c_str();
}
int GetParticles(uint16_t* data, int max_size) {
  int i = 0;
  for (auto p : g_scene.particles.p) {
    if (i + 1 >= max_size) {
      break;
    }
    data[i] = (1 + p.x) * gameinst->width_ * 0.5 / kScale;
    data[i + 1] = (1 - p.y) * gameinst->height_ * 0.5 / kScale;
    i += 2;
  }
  return i;
}
} // extern "C"

int main() {
  state_pause = false;
  state_quit = false;
  const int width = 800;
  const int height = 800;
  gameinst = std::unique_ptr<Game>(new Game(width, height));
  gameinst->partsys->RemoveLastPortal();
  g_view = std::make_unique<ViewText>();
  g_view->SetScene(g_scene);
  emscripten_set_canvas_element_size(
      "#canvas", width / kScale, height / kScale);
  emscripten_set_main_loop(main_loop, 30, 1);
}
