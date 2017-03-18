/*
  RENDERER
*/

#include "geometry.hpp"
#include "particles_system.hpp"
#include <GL/glut.h>
#include <iostream>
#include <functional>
class renderer
{
protected:
  particles_system* PS;
  projection PR;
public:
  renderer(particles_system* _PS, const projection& _PR) : PS(_PS), PR(_PR) {;}
  virtual void draw_particles() = 0;
  virtual void draw_frame() = 0;
};

class renderer_opengl : public renderer
{
public:
  void draw_circle(int x, int y, int r, rgb color)
  {
    glColor3f(color.r, color.g, color.b);

    float delta_theta = 0.1;

    //glBegin( GL_POLYGON );
    glBegin( GL_LINE_LOOP );

    for( float angle = 0; angle < 2*PI; angle += delta_theta )
    glVertex3f( r*cos(angle)+x, r*sin(angle)+y, 0 );

    glEnd();
  }
  void draw_circle(mindex c, int r, rgb color)
  {
    draw_circle(c.i, c.j, r, color);
  }
  void draw_line(mindex A, mindex B)
  {
    glColor3f(1., 1., 1.);
    glBegin( GL_LINES ); // OR GL_LINE_LOOP

    glVertex3f( A.i, A.j, 0 );
    glVertex3f( B.i, B.j, 0 );

    glEnd();
  }
  void draw_particles()
  {
    for(std::size_t k=0; k<PS->P.size(); ++k)
    {
      auto& part=PS->P[k];
      vect p=part.p;
      int r=PR.convert(vect(part.r,0.)).i-PR.convert(vect(0.,0.)).i;
      draw_circle(PR.convert(p), r, part.color);
    }
  }
  void draw_frame()
  {
    mindex A=PR.Rscreen.A, B=PR.Rscreen.B;
    draw_line(mindex(A.i,A.j), mindex(B.i,A.j));
    draw_line(mindex(A.i,B.j), mindex(B.i,B.j));
    draw_line(mindex(A.i,A.j), mindex(A.i,B.j));
    draw_line(mindex(B.i,A.j), mindex(B.i,B.j));
  }
  renderer_opengl(particles_system* _PS, const projection& _PR) : renderer(_PS, _PR)
  {

  }
};
