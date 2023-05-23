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

  Particles particles;
  span<std::array<Portal, 2>> portals;
};

class View {
 public:
  virtual ~View() = default;
  virtual const Scene& GetScene() const = 0;
  virtual void SetScene(const Scene&) = 0;
  virtual void Draw() = 0;
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
    std::cout << "particles.p\n";
    for (auto p : scene_.particles.p) {
      std::cout << p << ' ';
    }
    std::cout << '\n';

    std::cout << "particles.v\n";
    for (auto p : scene_.particles.v) {
      std::cout << p << ' ';
    }
    std::cout << '\n';

    std::cout << "portals\n";
    for (auto p : scene_.portals) {
      std::cout << p[0].pa << ' ' << p[0].pb << ' ' //
                << p[1].pa << ' ' << p[1].pb << '\n';
    }
    std::cout << '\n';
  }

 private:
  Scene scene_;
};
