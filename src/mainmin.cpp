#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "game.h"
#include "logger.h"
#include "scene.h"
#include "view_gl.h"

std::chrono::milliseconds last_frame_time;
Scal last_frame_game_time;
std::chrono::milliseconds last_report_time;
std::atomic<Scal> next_game_time_target;

std::unique_ptr<Game> gameinst;

std::atomic<bool> flag_display;
std::atomic<bool> state_pause;
std::atomic<bool> state_quit;

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
              << ", t=" << gameinst->partsys->GetTime()
              << ", steps=" << gameinst->partsys->GetNumSteps() << std::endl;
    last_report_time = new_frame_time;
  }

  const Scal game_rate_target = 1.5;
  if (!state_pause) {
    next_game_time_target = new_frame_game_time + game_rate_target / fps;
  }

  last_frame_time = new_frame_time;
  last_frame_game_time = new_frame_game_time;

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

  {
    std::lock_guard<std::mutex> lg(gameinst->partsys->m_buffer_);
    gameinst->partsys->SetRendererReadyForNext(true);
  }
  flag_display = false;
}

void cycle() {
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_nested(0);
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

  while (!state_quit) {
    gameinst->partsys->step(next_game_time_target, state_pause);
    if (state_pause) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  }
  std::cout << "Computation finished" << std::endl;
}

int main() {
  last_frame_game_time = 0.;
  state_pause = false;
  state_quit = false;

  int width = kInitWidth;
  int height = kInitHeight;

  gameinst = std::unique_ptr<Game>(new Game(width, height));

  // Create view.
  g_view = std::make_unique<ViewGl>(
      gameinst.get(), gameinst->partsys.get(), width, height, state_quit,
      state_pause);
  g_view->SetScene(g_scene);

  std::thread computation_thread(cycle);

  while (!state_quit) {
    g_view->Control();
    display();
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }

  computation_thread.join();
}
