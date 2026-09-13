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
#include <cstdlib>
#include <set>
#include <string>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif
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

  UpdateScene(*gameinst->partsys, g_data, g_scene);
  g_view->SetScene(g_scene);
  g_view->Draw();
  flag_display = false;
}

#ifdef _OPENMP
// Number of physical cores, or 0 if unknown.
static size_t GetNumCores() {
#if defined(__APPLE__)
  int value = 0;
  size_t size = sizeof(value);
  if (sysctlbyname("hw.physicalcpu", &value, &size, nullptr, 0) == 0 &&
      value > 0) {
    return static_cast<size_t>(value);
  }
  return 0;
#elif defined(__linux__)
  // Every processor entry reports the core it runs on, so the number of
  // distinct cores is the number of entries without the ones sharing a core.
  std::ifstream fin("/proc/cpuinfo");
  std::set<std::pair<int, int>> cores;
  int package = 0;
  std::string line;
  while (std::getline(fin, line)) {
    const size_t colon = line.find(':');
    if (colon == std::string::npos) {
      continue;
    }
    const std::string key = line.substr(0, colon);
    const int value = std::atoi(line.c_str() + colon + 1);
    if (key.compare(0, 11, "physical id") == 0) {
      package = value;
    } else if (key.compare(0, 7, "core id") == 0) {
      cores.insert({package, value});
    }
  }
  return cores.size();
#else
  return 0;
#endif
}

// OpenMP starts one thread per logical CPU, but the block loops in
// Particles::step() are memory-bound and gain nothing from a second thread on
// the same core, while the extra threads add synchronization at every barrier
// and take time from the renderer. Use one thread per core unless the
// environment asks for a number.
static void SetDefaultNumThreads() {
  if (!std::getenv("OMP_NUM_THREADS")) {
    const size_t cores = GetNumCores();
    if (cores > 0 && cores < static_cast<size_t>(omp_get_max_threads())) {
      omp_set_num_threads(static_cast<int>(cores));
    }
  }
  std::cout << "OpenMP threads: " << omp_get_max_threads() << std::endl;
}
#endif

int main() {
#ifdef _OPENMP
  SetDefaultNumThreads();
#endif
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
