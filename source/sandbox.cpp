#include <array>
#include <iostream>
#include <vector>

#include "scene.h"
#include "span.h"

struct Particles {
  std::vector<Vect> p;
  std::vector<Vect> v;
};

int main() {
  const size_t n = 3;
  Particles particles;
  for (size_t i = 0; i < n; ++i) {
    particles.p.push_back(i);
    particles.v.push_back(i * 0.1);
  }
  using Portal = Scene::Portal;
  std::vector<std::array<Portal, 2>> portals;
  for (Scal a : {0., 0.5}) {
    portals.push_back({
        Portal{Vect(a, 0), Vect(a, 1)},
        Portal{Vect(a + 1, 0), Vect(a + 1, 1)},
    });
  }

  Scene s;
  s.particles.p = particles.p;
  s.particles.v = particles.v;
  s.portals = portals;

  TextView v;
  v.SetScene(s);
  v.Draw();
}
