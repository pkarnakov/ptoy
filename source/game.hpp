/*
  GAME
*/

#include <vector>
#include "particles_system.hpp"
#include "renderer.hpp"
#include "geometry.hpp"
#include <mutex>

class game
{
public:
  std::unique_ptr<particles_system> PS;
  std::unique_ptr<renderer> R;
  game(int width, int height)
  {
    PS=std::unique_ptr<particles_system>(new particles_system);
    R=std::unique_ptr<renderer>(new renderer_opengl(PS.get()));
    SetWindowSize(width, height);
  }
  void SetWindowSize(int width, int height) {
    vect A(-1.,-1.), B(-1 + 2. * width / 800, -1. + 2. * height / 800);
    PS->SetDomain(rect_vect(A, B));
  }
};
