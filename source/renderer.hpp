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
public:
  renderer(particles_system* _PS) : PS(_PS) {;}
  virtual void draw_particles() = 0;
  virtual void draw_frame() = 0;
};

class renderer_opengl : public renderer
{
public:
  void draw_circle(Scal x, Scal y, Scal r, rgb color)
  {
    glColor3f(color.r, color.g, color.b);

    const size_t num_segments = 8;
    const Scal theta = 2 * PI / num_segments;
    const Scal c = cos(theta);
    const Scal s = sin(theta);
    Scal dx = r;
    Scal dy = 0.;

    glBegin( GL_POLYGON );
    //glBegin( GL_LINE_LOOP );
    //glDrawArrays();

    for (size_t i = 0; i < num_segments; ++i) {
      glVertex3f(x + dx, y + dy, 0);
      const Scal dx_tmp = dx;
      dx = c * dx - s * dy;
      dy = s * dx_tmp + c * dy;
    }

    glEnd();
  }
  void draw_circle(vect c, Scal r, rgb color)
  {
    draw_circle(c.x, c.y, r, color);
  }
  void draw_circle(mindex c, int r, rgb color)
  {
    draw_circle(c.i, c.j, r, color);
  }
  void draw_line(vect A, vect B)
  {
    glColor3f(1., 1., 1.);
    glBegin( GL_LINES ); // OR GL_LINE_LOOP

    glVertex2f(A.x, A.y);
    glVertex2f(B.x, B.y);

    glEnd();
  }
  void draw_line(mindex A, mindex B)
  {
    glColor3f(1., 1., 1.);
    glBegin( GL_LINES ); // OR GL_LINE_LOOP

    glVertex2f(A.i, A.j);
    glVertex2f(B.i, B.j);

    glEnd();
  }
  void draw_particles()
  {
    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    rect_vect R = PS->GetDomain();
    glOrtho(R.A.x, R.B.x, R.A.y, R.B.y, -1.f, 1.f);
    auto particles = PS->GetParticles();
    for(std::size_t k=0; k<particles.size(); ++k)
    {
      auto& part=particles[k];
      //vect p=part.p;
      Scal f = 0.5 + part.v.length() / 7.; // color intensity
      f = std::min<Scal>(std::max<Scal>(f, 0.), 1.);
      //auto c = part.color;
      //rgb color(c.r * f, c.g * f, c.b * f);
      rgb color(f, 0., 0.);
      //draw_circle(part.p, part.r, color);
      draw_circle(part.p, kRadius * 0.9, color);
    }
    glPopMatrix();
  }
  void draw_frame()
  {
    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    rect_vect R = PS->GetDomain();
    glLineWidth(5.0);
    vect A = R.A;
    vect B = R.B;
    draw_line(vect(A.x,A.y), vect(B.x,A.y));
    draw_line(vect(A.x,B.y), vect(B.x,B.y));
    draw_line(vect(A.x,A.y), vect(A.x,B.y));
    draw_line(vect(B.x,A.y), vect(B.x,B.y));
    glPopMatrix();
  }
  renderer_opengl(particles_system* _PS) : renderer(_PS)
  {

  }
};
