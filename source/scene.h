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

  Particles particles;
  span<std::array<Portal, 2>> portals;
};

class View {
 public:
  virtual const Scene& GetScene() const = 0;
  virtual void SetScene(const Scene&) = 0;
  virtual void Draw(const Scene&) = 0;
};

class TextView : public View {
 public:
  const Scene& GetScene() const override {
    return scene_;
  }
  void SetScene(const Scene& scene) override {
    scene_ = scene;
  }
  void Draw() override {
    std::cout << scene_.particles.p.size() << std::endl;
  }

 private:
  Scene scene_;
}
