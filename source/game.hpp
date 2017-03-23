/*
  GAME
*/

#include <vector>
#include "particles_system.hpp"
#include "renderer.hpp"
#include "geometry.hpp"

class game
{
public:
  std::unique_ptr<particles_system> PS;
  std::unique_ptr<renderer> R;
  game()
  {
    PS=std::unique_ptr<particles_system>(new particles_system);

    vect A(-1.,-1.), B(1.,1.);

    R=std::unique_ptr<renderer>(new renderer_opengl(PS.get()));

    double eps=0.01;
    // place env_objects
    PS->AddEnvObj(new line(vect(A.x, A.y), vect(B.x, A.y), eps));
    PS->AddEnvObj(new line(vect(A.x, B.y), vect(B.x, B.y), eps));
    PS->AddEnvObj(new line(vect(A.x, A.y), vect(A.x, B.y), eps));
    PS->AddEnvObj(new line(vect(B.x, A.y), vect(B.x, B.y), eps));
  }
};
