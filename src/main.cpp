#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "game.h"
#include "logger.h"
#include "macros.h"
#include "scene.h"

#if USEFLAG(BACKEND_TEXT)
#include "view_text.h"
#endif

#if USEFLAG(BACKEND_SDL)
#include "view_gl.h"
#endif

std::chrono::milliseconds g_last_wtime{0};
std::chrono::milliseconds g_report_wtime;
Scal g_last_gtime;

std::unique_ptr<Game> gameinst;

bool flag_display;
bool state_pause;
bool state_quit;

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

Scal Clip(Scal a, Scal lower, Scal upper) {
  return std::max(lower, std::min(upper, a));
}

void display() {
  if (flag_display) return;
  flag_display = true;

  const auto curr_wtime = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now().time_since_epoch());
  const Scal curr_gtime = gameinst->partsys->GetTime();

  if (!g_last_wtime.count()) {
    std::cout << "first frame\n";
    g_last_wtime = curr_wtime;
    g_report_wtime = curr_wtime;
  }

  const Scal frame_wtime = (curr_wtime - g_last_wtime).count() / 1000.;

  if ((curr_wtime - g_report_wtime).count() > 1000) {
    std::cout << "fps=" << 1 / frame_wtime
              << "  particles=" << gameinst->partsys->GetNumParticles()
              << "  speed=" << std::setprecision(3)
              << (curr_gtime - g_last_gtime) / frame_wtime
              << "  t=" << gameinst->partsys->GetTime() << std::endl;
    g_report_wtime = curr_wtime;
  }

  const Scal speed_target = 1.5;
  if (!state_pause) {
    const Scal next_game_time_target =
        curr_gtime + Clip(speed_target * frame_wtime, 0.02, 0.04);
    gameinst->partsys->SetRendererReadyForNext(true);
    gameinst->partsys->step(next_game_time_target, state_pause);
  }

  g_last_wtime = curr_wtime;
  g_last_gtime = curr_gtime;

  // Update scene.
  {
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
  }
  g_view->Draw();
  flag_display = false;
}

int main() {
  g_last_gtime = 0;
  state_pause = false;
  state_quit = false;

  int width = kInitWidth;
  int height = kInitHeight;

  gameinst = std::unique_ptr<Game>(new Game(width, height));

#if USEFLAG(BACKEND_TEXT)
  g_view = std::make_unique<ViewText>();
#endif

#if USEFLAG(BACKEND_SDL)
  g_view = std::make_unique<ViewGl>(
      gameinst.get(), gameinst->partsys.get(), width, height, state_quit,
      state_pause);
#endif

  g_view->SetScene(g_scene);

  while (!state_quit) {
    g_view->Control();
    display();
  }
}
