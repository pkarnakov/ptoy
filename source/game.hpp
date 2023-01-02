#include <mutex>
#include <vector>
#include "geometry.hpp"
#include "particles_system.hpp"
#include "renderer.hpp"

constexpr int kInitWidth = 800;
constexpr int kInitHeight = 800;

class Game {
 public:
  std::unique_ptr<particles_system> partsys;
  std::unique_ptr<renderer> rendinst;
  int width_;
  int height_;
  Game(int width, int height) {
    partsys = std::unique_ptr<particles_system>(new particles_system);
    rendinst = std::unique_ptr<renderer>(new renderer(partsys.get()));
    SetWindowSize(width, height);
  }
  void SetWindowSize(int width, int height) {
    width_ = width;
    height_ = height;
    vect A(-1., -1.),
        B(-1 + 2. * width / kInitWidth, -1. + 2. * height / kInitHeight);
    partsys->PushResize(rect_vect(A, B));
    rendinst->SetWindowSize(width, height);
  }
};
