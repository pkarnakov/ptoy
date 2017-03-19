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

    const size_t num_segments = 8;
    const float theta = 2 * PI / num_segments;
    const float c = cos(theta);
    const float s = sin(theta);
    float dx = r;
    float dy = 0.;

    glBegin( GL_POLYGON );
    //glBegin( GL_LINE_LOOP );
    //glDrawArrays();

    for (size_t i = 0; i < num_segments; ++i) {
      glVertex3f(x + dx, y + dy, 0);
      const float dx_tmp = dx;
      dx = c * dx - s * dy;
      dy = s * dx_tmp + c * dy;
    }

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
    auto particles = PS->GetParticles();
    for(std::size_t k=0; k<particles.size(); ++k)
    {
      auto& part=particles[k];
      vect p=part.p;
      double f = 0.5 + part.v.length() / 7.; // color intensity
      f = std::min(1., std::max(0., f));
      auto c = part.color;
      rgb color(c.r * f, c.g * f, c.b * f);
      int r = PR.convert(vect(part.r,0.)).i-PR.convert(vect(0.,0.)).i;
      draw_circle(PR.convert(p), r, color);
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
