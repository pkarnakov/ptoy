#include <SDL_opengl.h>
#include <functional>
#include <iostream>

#include "geometry.hpp"
#include "particles_system.hpp"
#include "logger.h"

#define CHECK_ERROR()                                         \
  do {                                                        \
    GLenum error = glGetError();                              \
    fassert(                                                  \
        error == GL_NO_ERROR,                                 \
        reinterpret_cast<const char*>(gluErrorString(error))) \
  } while (0);

class renderer {
  int width_, height_;
  particles_system* PS;

 public:
  void SetWindowSize(int width, int height) {
    width_ = width;
    height_ = height;
  }
  void draw_circle(Scal x, Scal y, Scal r, rgb color) {
    glColor3f(color.r, color.g, color.b);

    const size_t num_segments = 8;
    const Scal theta = 2 * PI / num_segments;
    const Scal c = cos(theta);
    const Scal s = sin(theta);
    Scal dx = r;
    Scal dy = 0.;

    glBegin(GL_POLYGON);

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
  void draw_line(vect A, vect B, rgb color, Scal width) {
    glColor3f(color.r, color.g, color.b);
    glBegin(GL_POLYGON);

    const vect n = vect(A.y - B.y, B.x - A.x).GetNormalized();
    auto v = [](vect p) { //
      glVertex2f(p.x, p.y);
    };
    v(A + n * (width * 0.5));
    v(A - n * (width * 0.5));
    v(B - n * (width * 0.5));
    v(B + n * (width * 0.5));

    glEnd();
  }
  void draw_line(vect A, vect B, Scal width) {
    draw_line(A, B, rgb(1., 1., 1.), width);
  }
  void draw_particles() {
    auto particles = PS->GetParticles();
    for (std::size_t k = 0; k < particles.size(); ++k) {
      auto& part = particles[k];
      // vect p=part.p;
      Scal f = 0.5 + part.v.length() / 7.; // color intensity
      f = std::min<Scal>(std::max<Scal>(f, 0.), 1.);
      // auto c = part.color;
      // rgb color(c.r * f, c.g * f, c.b * f);
      rgb color(f, 0., 0.);
      // draw_circle(part.p, part.r, color);
      draw_circle(part.p, kRadius * 0.9, color);
    }
  }
  void draw_frame() {
    rect_vect R = PS->GetDomain();
    vect dA = R.A;
    vect dB = R.B;
    const Scal width = 0.005;
    draw_line(vect(dA.x, dA.y), vect(dB.x, dA.y), width);
    draw_line(vect(dA.x, dB.y), vect(dB.x, dB.y), width);
    draw_line(vect(dA.x, dA.y), vect(dA.x, dB.y), width);
    draw_line(vect(dB.x, dA.y), vect(dB.x, dB.y), width);
  }
  void DrawBonds() {
    const auto& nr = PS->GetNoRendering();
    const auto& data = PS->GetBlockData();
    const auto& bbi = PS->GetBlockById();
    for (auto bond : PS->GetBonds()) {
      const auto& a = bbi[bond.first];
      const auto& b = bbi[bond.second];
      assert(a.first != blocks::kBlockNone);
      assert(b.first != blocks::kBlockNone);
      if (!nr.count(bond)) {
        const Scal width = 0.0075;
        draw_line(
            data.position[a.first][a.second], data.position[b.first][b.second],
            width);
      }
    }
  }
  void DrawFrozen() {
    const auto& data = PS->GetBlockData();
    const auto& bbi = PS->GetBlockById();
    for (auto id : PS->GetFrozen()) {
      const auto& a = bbi[id];
      assert(a.first != blocks::kBlockNone);
      draw_circle(
          data.position[a.first][a.second], kRadius * 0.5, rgb(0., 0., 0.));
    }
  }
  void DrawPortals() {
    const auto& portals = PS->GetPortals();
    const Scal width = 0.0075;
    for (auto& pair : portals) {
      draw_line(pair[0].begin, pair[0].end, rgb(0., .5, 1.), width);
      draw_line(pair[1].begin, pair[1].end, rgb(1., .5, 0.), width);
    }
    if (PS->portal_stage_ == 0) {
      if (PS->portal_mouse_moving_) {
        draw_line(
            PS->portal_begin_, PS->portal_current_, rgb(0., .5, 1.), width);
      }
    } else {
      draw_line(
          PS->portal_prev_.first, PS->portal_prev_.second, rgb(0., .5, 1.),
          width);
      if (PS->portal_mouse_moving_) {
        draw_line(
            PS->portal_begin_, PS->portal_current_, rgb(1., .5, 0.), width);
      }
    }
  }
  void DrawAll() {
    //draw_particles(); // done by shaders in main.cpp
    draw_frame();
    DrawBonds();
    DrawFrozen();
    DrawPortals();
  }
  renderer(particles_system* _PS) : PS(_PS) {}
};
