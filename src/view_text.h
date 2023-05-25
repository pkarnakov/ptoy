#pragma once

#include <iostream>

#include "scene.h"

class ViewText : public View {
 public:
  ViewText() = default;
  const Scene& GetScene() const override {
    return scene_;
  }
  void SetScene(const Scene& scene) override {
    scene_ = scene;
  }
  void Control() override {}
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
