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
  int width_, height_;
  particles_system* PS;
public:
  void SetWindowSize(int width, int height) {
    width_ = width;
    height_ = height;
  }
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
  void draw_circle(vect c, Scal r, rgb color) {
    draw_circle(c.x, c.y, r, color);
  }
  void draw_circle(mindex c, int r, rgb color) {
    draw_circle(c.i, c.j, r, color);
  }
  void draw_line(vect A, vect B, rgb color) {
    glColor3f(color.r, color.g, color.b);
    glBegin(GL_LINES); // OR GL_LINE_LOOP

    glVertex2f(A.x, A.y);
    glVertex2f(B.x, B.y);

    glEnd();
  }
  void draw_line(vect A, vect B) {
    draw_line(A, B, rgb(1., 1., 1.));
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
    vect A(-1.,-1.), B(-1 + 2. * width_ / 800, -1. + 2. * height_ / 800);
    glOrtho(A.x, B.x, A.y, B.y, -1.f, 1.f);
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
    vect A(-1.,-1.), B(-1 + 2. * width_ / 800, -1. + 2. * height_ / 800);
    glOrtho(A.x, B.x, A.y, B.y, -1.f, 1.f);
    rect_vect R = PS->GetDomain();
    glLineWidth(5.0);
    vect dA = R.A;
    vect dB = R.B;
    draw_line(vect(dA.x,dA.y), vect(dB.x,dA.y));
    draw_line(vect(dA.x,dB.y), vect(dB.x,dB.y));
    draw_line(vect(dA.x,dA.y), vect(dA.x,dB.y));
    draw_line(vect(dB.x,dA.y), vect(dB.x,dB.y));
    glPopMatrix();
  }
  void DrawBonds() {
    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    vect A(-1.,-1.), B(-1 + 2. * width_ / 800, -1. + 2. * height_ / 800);
    glOrtho(A.x, B.x, A.y, B.y, -1.f, 1.f);
    glLineWidth(3.0);
    const auto& data = PS->GetBlockData();
    const auto& bbi = PS->GetBlockById();
    for (auto bond : PS->GetBonds()) {
      const auto& a = bbi[bond.first];
      const auto& b = bbi[bond.second];
      assert(a.first != blocks::kBlockNone);
      assert(b.first != blocks::kBlockNone);
      draw_line(
          data.position[a.first][a.second],
          data.position[b.first][b.second]);
    }
    glPopMatrix();
  }
  void DrawFrozen() {
    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    vect A(-1.,-1.), B(-1 + 2. * width_ / 800, -1. + 2. * height_ / 800);
    glOrtho(A.x, B.x, A.y, B.y, -1.f, 1.f);
    glLineWidth(3.0);
    const auto& data = PS->GetBlockData();
    const auto& bbi = PS->GetBlockById();
    for (auto id : PS->GetFrozen()) {
      const auto& a = bbi[id];
      assert(a.first != blocks::kBlockNone);
      draw_circle(data.position[a.first][a.second], 
                 kRadius * 0.5, rgb(0.,0.,0.));
    }
    glPopMatrix();
  }
  void DrawPortals() {
    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    vect A(-1.,-1.), B(-1 + 2. * width_ / 800, -1. + 2. * height_ / 800);
    glOrtho(A.x, B.x, A.y, B.y, -1.f, 1.f);
    glLineWidth(3.0);
    const auto& portals = PS->GetPortals();
    for (auto& pair : portals) {
      draw_line(pair.first.begin, pair.first.end, rgb(0., .5, 1.));
      draw_line(pair.second.begin, pair.second.end, rgb(1., .5, 0.));
    }
    if (PS->portal_stage_ == 0) {
      if (PS->portal_mouse_moving_) {
        draw_line(PS->portal_begin_, PS->portal_current_, rgb(0., .5, 1.));
      }
    } else  {
      draw_line(PS->portal_prev_.first, 
                PS->portal_prev_.second, rgb(0., .5, 1.));
      if (PS->portal_mouse_moving_) {
        draw_line(PS->portal_begin_, PS->portal_current_, rgb(1., .5, 0.));

      }
    }
    glPopMatrix();
  }
  void DrawAll() {
    draw_particles();
    draw_frame();
    DrawBonds();
    DrawFrozen();
    DrawPortals();
  }
  renderer(particles_system* _PS) : PS(_PS) {}
};
