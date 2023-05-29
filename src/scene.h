#pragma once

#include "geometry.h"
#include "span.h"

struct Scene {
  struct Particles {
    span<Vect> p; // Position.
    span<Vect> v; // Velocity.
  };
  struct Portal {
    Vect pa;
    Vect pb;
  };
  struct Bond {
    Vect pa;
    Vect pb;
  };

  Particles particles;
  span<std::array<Portal, 2>> portals;
  span<Bond> bonds;
  span<Vect> frozen;
};

class View {
 public:
  virtual ~View() = default;
  virtual const Scene& GetScene() const = 0;
  virtual void SetScene(const Scene&) = 0;
  virtual void Control() = 0;
  virtual void Draw() = 0;
};
