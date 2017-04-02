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
  std::unique_ptr<renderer_opengl> R;
  int width_;
  int height_;
  game(int width, int height)
  {
    PS=std::unique_ptr<particles_system>(new particles_system);
    R=std::unique_ptr<renderer_opengl>(new renderer_opengl(PS.get()));
    SetWindowSize(width, height);
  }
  void SetWindowSize(int width, int height) {
    width_ = width;
    height_ = height;
    vect A(-1.,-1.), B(-1 + 2. * width / 800, -1. + 2. * height / 800);
    PS->PushResize(rect_vect(A, B));
    R->SetWindowSize(width, height);
  }
};
