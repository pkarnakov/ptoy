#include <mutex>
#include <vector>
#include "geometry.h"
#include "particles.h"

constexpr int kInitWidth = 800;
constexpr int kInitHeight = 800;

class Game {
 public:
  std::unique_ptr<Particles> partsys;
  int width_;
  int height_;
  Game(int width, int height) {
    partsys = std::unique_ptr<Particles>(new Particles);
    SetWindowSize(width, height);
  }
  void SetWindowSize(int width, int height) {
    width_ = width;
    height_ = height;
    Vect A(-1., -1.),
        B(-1 + 2. * width / kInitWidth, -1. + 2. * height / kInitHeight);
    partsys->PushResize(RectVect(A, B));
  }
};
