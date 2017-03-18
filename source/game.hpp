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

    R=std::unique_ptr<renderer>(new renderer_opengl(PS.get(), projection( rect_vect(A, B), 
                                                                    rect_mindex(mindex(000,000), mindex(500,500)) ) ));

    double eps=0.01;
    // place env_objects
    PS->ENVOBJ.push_back(std::unique_ptr<env_object>(new line(vect(A.x, A.y), vect(B.x, A.y), eps)));
    PS->ENVOBJ.push_back(std::unique_ptr<env_object>(new line(vect(A.x, B.y), vect(B.x, B.y), eps)));
    PS->ENVOBJ.push_back(std::unique_ptr<env_object>(new line(vect(A.x, A.y), vect(A.x, B.y), eps)));
    PS->ENVOBJ.push_back(std::unique_ptr<env_object>(new line(vect(B.x, A.y), vect(B.x, B.y), eps)));
  }
};