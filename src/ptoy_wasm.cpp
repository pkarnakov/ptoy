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

const int g_width = 800;
const int g_height = 800;
Scene g_scene;
SceneData g_data;
std::shared_ptr<Game> g_gameinst;
std::shared_ptr<Control> g_control;
std::string g_buf;
bool state_pause;

static void main_loop() {
  auto gameinst = g_gameinst;
  gameinst->partsys->SetRendererReadyForNext(true);
  const auto dt = 0.03;
  gameinst->partsys->step(gameinst->partsys->GetTime() + dt, state_pause);
  UpdateScene(*gameinst->partsys, g_data, g_scene);
  EM_ASM_({ draw(); });
}

extern "C" {
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
int GetBonds(uint16_t* data, int max_size) {
  auto gameinst = g_gameinst;
  const int entrysize = 4;
  int i = 0;
  for (auto bond : g_scene.bonds) {
    if (i + entrysize > max_size) {
      break;
    }
    auto append = [&](Vect p) {
      data[i + 0] = (1 + p.x) * gameinst->width_ * 0.5;
      data[i + 1] = (1 - p.y) * gameinst->height_ * 0.5;
      i += 2;
    };
    append(bond.pa);
    append(bond.pb);
  }
  return i;
}
int GetFrozen(uint16_t* data, int max_size) {
  auto gameinst = g_gameinst;
  const int entrysize = 2;
  int i = 0;
  for (auto p : g_scene.frozen) {
    if (i + entrysize > max_size) {
      break;
    }
    auto append = [&](Vect p) {
      data[i + 0] = (1 + p.x) * gameinst->width_ * 0.5;
      data[i + 1] = (1 - p.y) * gameinst->height_ * 0.5;
      i += 2;
    };
    append(p);
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
  UpdateScene(*g_gameinst->partsys, g_data, g_scene);
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
void SetGravityVect(float gx, float gy) {
  auto gameinst = g_gameinst;
  gameinst->partsys->SetGravityVect({gx, gy});
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
  emscripten_set_main_loop(main_loop, 0, 1);
}
