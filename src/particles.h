#pragma once

#include <functional>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "blocks.h"
#include "geometry.h"

extern const Scal kRadius;
extern const Scal kPortalThickness;

struct Particle {
  Vect p;
  Vect v;
  Particle() {}
  Particle(Vect p_, Vect v_) : p(p_), v(v_) {}
};

Vect F12(Vect p1, Vect p2, Scal sigma, Scal R);
Vect F12wall(Vect p1, Vect p2);
Vect F12(Vect p1, Vect p2);

class Obstacle {
 public:
  virtual ~Obstacle() = default;
  virtual Vect F(Vect p, Vect v) = 0;
  virtual bool IsClose(Vect p, Scal R) = 0;
};

class Line : public Obstacle {
  Vect A, B;
  Vect GetNearest(Vect p) {
    Vect Q;
    Scal lambda = (B - A).dot(p - A) / (B - A).dot(B - A);
    if (lambda > 0. && lambda < 1.) {
      Q = A + (B - A) * lambda;
    } else {
      Q = (p.dist(A) < p.dist(B)) ? A : B;
    }
    return Q;
  }

 public:
  Line(Vect _A, Vect _B) : A(_A), B(_B) {}
  Vect F(Vect p, Vect) override {
    return F12wall(p, GetNearest(p));
  }
  bool IsClose(Vect p, Scal R) override {
    return GetNearest(p).dist(p) < R + kRadius;
  }
};

class Particles {
 public:
  Particles();
  ~Particles();
  struct Portal {
    Vect begin, end;
    std::vector<size_t> blocks;
    Vect GetNearest(Vect p) {
      const Vect A = begin, B = end;
      Vect Q;
      Scal lambda = (B - A).dot(p - A) / (B - A).dot(B - A);
      if (lambda > 0. && lambda < 1.) {
        Q = A + (B - A) * lambda;
      } else {
        Q = (p.dist(A) < p.dist(B)) ? A : B;
      }
      return Q;
    }
    bool IsClose(Vect p, Scal R) {
      return GetNearest(p).dist(p) < R + 2 * kPortalThickness;
    }
  };
  const std::vector<std::array<Portal, 2>>& GetPortals() const {
    return portals_;
  }
  void ApplyPortals();
  void ApplyPortalsForces();
  void DetectPortals();
  void MoveToPortal(
      Vect& position, Vect& velocity, const Portal& src, const Portal& dest);
  void SetParticleBuffer();
  void RemoveLastPortal() {
    remove_last_portal_ = true;
  }
  void AddEnvObj(Obstacle* env);
  void ClearEnvObj() {
    obstacles_.clear();
  }
  void UpdateEnvObj();
  void ResetEnvObjFrame(RectVect new_domain) {
    const Vect A = new_domain.A, B = new_domain.B;
    ClearEnvObj();
    AddEnvObj(new Line(Vect(A.x, A.y), Vect(B.x, A.y)));
    AddEnvObj(new Line(Vect(A.x, B.y), Vect(B.x, B.y)));
    AddEnvObj(new Line(Vect(A.x, A.y), Vect(A.x, B.y)));
    AddEnvObj(new Line(Vect(B.x, A.y), Vect(B.x, B.y)));
    UpdateEnvObj();
  }
  void SetDomain(RectVect new_domain) {
    domain = new_domain;
    blocks_.SetDomain(domain);
  }
  void PushResize(RectVect new_domain) {
    resize_queue_ = new_domain;
  }
  const std::set<std::pair<int, int>>& GetBonds() const {
    return bonds_;
  }
  const std::set<int>& GetFrozen() const {
    return frozen_;
  }
  void step(Scal time_target, bool quit);
  void step_old(Scal time_target, bool quit);
  void SetForce(Vect center, bool enabled);
  void SetForce(Vect center);
  void SetForce(bool enabled);
  void SetForceAttractive(bool value) {
    force_attractive_ = value;
  }
  void BondsStart(Vect point);
  void BondsMove(Vect point);
  void BondsStop(Vect point);
  void CheckBonds();
  void FreezeStart(Vect point);
  void FreezeMove(Vect point);
  void FreezeStop(Vect point);
  void PickStart(Vect point);
  void PickMove(Vect point);
  void PickStop(Vect point);
  void PortalStart(Vect point);
  void PortalMove(Vect point);
  void PortalStop(Vect point);
  Scal GetTime() const {
    return time_;
  }
  RectVect GetDomain() const {
    return domain;
  }

  // These are only for external use (TODO: check or ensure)
  bool GetGravity() const {
    return gravity_enable_;
  }
  void SetGravity(bool flag) {
    gravity_enable_ = flag;
  }
  Vect GetGravityVect() const {
    return gravity_;
  }
  void SetGravityVect(Vect gravity) {
    gravity_ = gravity;
  }
  const std::vector<Particle>& GetParticles() const {
    return particle_buffer_;
  }
  size_t GetNumParticles() const {
    return blocks_buffer_.GetNumParticles();
  }
  size_t GetNumPerCell() const {
    return blocks_buffer_.GetNumPerCell();
  }
  const std::vector<std::pair<size_t, size_t>> GetBlockById() const {
    return blocks_buffer_.GetBlockById();
  }
  const Blocks::BlockData& GetBlockData() const {
    return blocks_buffer_.GetData();
  };
  void SetRendererReadyForNext(bool value) {
    renderer_ready_for_next_ = value;
  }

 private:
  RectVect domain;
  RectVect resize_queue_;
  Blocks blocks_;
  Scal time_;
  Vect gravity_;
  void CalcForces(size_t i);
  void AppendBondForces();
  void ApplyFrozen();
  std::vector<std::unique_ptr<Obstacle>> obstacles_;
  std::vector<std::vector<size_t>> obstacles_in_block_;
  Vect force_center;
  bool force_enabled;
  bool force_attractive_ = false;
  bool gravity_enable_ = true;
  std::vector<Particle> particle_buffer_;
  Blocks blocks_buffer_;
  int bonds_prev_particle_id_;
  bool bonds_enabled_ = false;
  int pick_particle_id_;
  bool pick_enabled_ = false;
  Vect pick_pointer_;
  std::set<std::pair<int, int>> bonds_; // Id's of particles with bonds.
  std::set<int> frozen_; // Id's of frozen particles.
  bool freeze_enabled_ = false;
  int freeze_last_id_;
  bool renderer_ready_for_next_ = true;
  bool portal_enabled_ = false;
  bool remove_last_portal_;

 public:
  int portal_stage_ = 0;
  bool portal_mouse_moving_ = false;
  Vect portal_begin_;
  Vect portal_current_;
  std::pair<Vect, Vect> portal_prev_;
  std::vector<std::array<Portal, 2>> portals_;
  std::vector<int> particle_to_move_; // 1: to move, 0: otherwise
  void UpdatePortalBlocks(Portal& portal);
  std::set<std::pair<int, int>> no_rendering_;
  std::set<std::pair<int, int>> no_rendering_buffer_;
  const std::set<std::pair<int, int>>& GetNoRendering() const {
    return no_rendering_buffer_;
  }
};
