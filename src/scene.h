#pragma once

#include <vector>

#include "geometry.h"
#include "span.h"

class Particles;

// Everything a view needs to draw one frame. Holds no data of its own, only
// spans into a SceneData filled by UpdateScene().
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
  // Pairs of portals. A portal with pa == pb is absent, which happens for the
  // pair currently being drawn.
  span<std::array<Portal, 2>> portals;
  span<Bond> bonds;
  span<Vect> frozen;
  RectVect domain;
};

// Storage for the data referenced by Scene.
struct SceneData {
  std::vector<Vect> position;
  std::vector<Vect> velocity;
  std::vector<std::array<Scene::Portal, 2>> portals;
  std::vector<Scene::Bond> bonds;
  std::vector<Vect> frozen;
};

// Collects the current state of `partsys` in `data` and points `scene` to it.
void UpdateScene(const Particles& partsys, SceneData& data, Scene& scene);

class View {
 public:
  virtual ~View() = default;
  virtual const Scene& GetScene() const = 0;
  virtual void SetScene(const Scene&) = 0;
  virtual void Control() = 0;
  virtual void Draw() = 0;
};
