#pragma once

#include "game.h"
#include "scene.h"

class ViewGl : public View {
 public:
  ViewGl(
      Game* gameinst_, Particles* partsys_, unsigned width_, unsigned height_,
      std::atomic<bool>& state_quit_, std::atomic<bool>& state_pause_);
  ~ViewGl();
  const Scene& GetScene() const override {
    return scene_;
  }
  void SetScene(const Scene& scene) override {
    scene_ = scene;
  }
  void Control() override;
  void Draw() override;

 private:
  Scene scene_;

  struct Imp;
  std::unique_ptr<Imp> imp;
};
