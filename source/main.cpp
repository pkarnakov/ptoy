#include <iostream>

#include <SDL.h>
#include "game.hpp"

#include <atomic>
#include <chrono>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

// TODO: reorder header includes (std first)
// TODO: fix formatting

using namespace std::chrono;

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

  // if(time_past_from_last_frame<milliseconds(100)) return;

  milliseconds frame_duration(int(1000.0 / fps));
  milliseconds time_residual = frame_duration - time_past_from_last_frame;

  // cout<<"sleep for "<<time_residual.count()<<endl;
  // std::this_thread::sleep_for(time_residual);
  // std::this_thread::sleep_for(frame_duration);

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

  const Scal game_rate_target = 10.;
  // const Scal game_rate_target = 1.;
  if (!pause) {
    next_game_time_target = new_frame_game_time + game_rate_target / fps;
  }

  last_frame_time = new_frame_time;
  last_frame_game_time = new_frame_game_time;
  // cout<<"Frame "<<frame_number++<<endl;

  /* clear the screen to white */
  glClear(GL_COLOR_BUFFER_BIT);

  {
    std::lock_guard<std::mutex> lg(G->PS->m_buffer_);
    G->R->DrawAll();
    G->PS->SetRendererReadyForNext(true);
  }

  glFlush();

  // glutSwapBuffers();

  flag_display = false;
}

void cycle() {
  // omp_set_dynamic(0);
  // omp_set_nested(0);
  // omp_set_num_threads(std::thread::hardware_concurrency());
#ifdef _OPENMP
  omp_set_num_threads(2);
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
  glClearColor(gray, gray, gray, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  SDL_GL_SwapWindow(window);

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
