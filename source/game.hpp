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
    rect_vect R(A, B);
    PS->SetDomain(R);

    std::lock_guard<std::mutex> lg(PS->m_ENVOBJ);

    // place env_objects
    PS->ClearEnvObj();
    PS->AddEnvObj(new line(vect(A.x, A.y), vect(B.x, A.y)));
    PS->AddEnvObj(new line(vect(A.x, B.y), vect(B.x, B.y)));
    PS->AddEnvObj(new line(vect(A.x, A.y), vect(A.x, B.y)));
    PS->AddEnvObj(new line(vect(B.x, A.y), vect(B.x, B.y)));
    PS->UpdateEnvObj();

  }
};
