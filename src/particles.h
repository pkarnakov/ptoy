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

// Presets for the grid of particles, in particles per side.
constexpr size_t kGridSmall = 10;
constexpr size_t kGridMedium = 25;
constexpr size_t kGridLarge = 45;
// Block added by AddParticleBlock(), in particles per side.
constexpr size_t kBlockGrid = 10;
// Where that block appears, relative to the domain.
constexpr Scal kBlockPositionX = 0.5;
constexpr Scal kBlockPositionY = 0.8;

struct particle {
  Vect p;
  Vect v;
  Vect f;
  particle() {}
  particle(Vect p_, Vect v_) : p(p_), v(v_) {}
};

Vect F12(Vect p1, Vect p2, Scal sigma, Scal R);
Vect F12wall(Vect p1, Vect p2);
Vect F12(Vect p1, Vect p2);

class env_object {
 public:
  virtual ~env_object() = default;
  virtual Vect F(Vect p, Vect v) = 0;
  virtual bool IsClose(Vect p, Scal R) = 0;
};

class line : public env_object {
  Vect A, B;
  // Scal eps;
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
  line(Vect _A, Vect _B) : A(_A), B(_B) {
    ;
  }
  Vect F(Vect p, Vect /*v*/) override {
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
  void AddEnvObj(env_object* env);
  void ClearEnvObj() {
    ENVOBJ.clear();
  }
  void UpdateEnvObj();
  void ResetEnvObjFrame(RectVect new_domain) {
    const Vect A = new_domain.A, B = new_domain.B;
    ClearEnvObj();
    AddEnvObj(new line(Vect(A.x, A.y), Vect(B.x, A.y)));
    AddEnvObj(new line(Vect(A.x, B.y), Vect(B.x, B.y)));
    AddEnvObj(new line(Vect(A.x, A.y), Vect(A.x, B.y)));
    AddEnvObj(new line(Vect(B.x, A.y), Vect(B.x, B.y)));
    UpdateEnvObj();
  }
  void SetDomain(RectVect new_domain) {
    domain = new_domain;
    Blocks.SetDomain(domain);
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
  // Replaces all particles with a grid of `size` by `size` of them.
  // Resets the simulation, see the definition.
  void SetParticleGrid(size_t size);
  // Adds a block of `size` by `size` particles centered at `relative`
  // position in the domain, where (0,0) is its lower left corner and (1,1) the
  // upper right one. Keeps the particles that are already there.
  void AddParticleBlock(Vect relative, size_t size);
  // Size of the grid passed to the last SetParticleGrid().
  size_t GetParticleGrid() const {
    return grid_size_;
  }
  void step(Scal time_target, bool quit);
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
  void CheckFrozen();
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
    return t;
  }
  size_t GetNumSteps() const {
    return static_cast<size_t>(t / dt);
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
  // Changes gravity by `steps`, negative pointing it further down and
  // positive further up, up to twice its default magnitude either way. Turns
  // gravity on, since it would have no effect otherwise.
  void ChangeGravity(int steps);
  void SetGravityVect(Vect gravity) {
    gravity_ = gravity;
  }
  const std::vector<particle>& GetParticles() const {
    return particle_buffer_;
  }
  size_t GetNumParticles() const {
    return blocks_buffer_.GetNumParticles();
  }
  size_t GetNumPerCell() const {
    return blocks_buffer_.GetNumPerCell();
  }
  const std::vector<std::pair<size_t, size_t>>& GetBlockById() const {
    return blocks_buffer_.GetBlockById();
  }
  const blocks::BlockData& GetBlockData() const {
    return blocks_buffer_.GetData();
  };
  void SetRendererReadyForNext(bool value) {
    renderer_ready_for_next_ = value;
  }
  // Portal pair currently being drawn with the mouse.
  struct PortalDrawing {
    int stage; // 0: drawing the first portal of a pair, 1: the second one.
    bool mouse_moving; // Whether a portal is being drawn right now.
    Vect begin, current; // Endpoints of the portal being drawn.
    std::pair<Vect, Vect> prev; // First portal of the pair, if stage is 1.
  };
  PortalDrawing GetPortalDrawing() const {
    return {portal_stage_, portal_mouse_moving_, portal_begin_, portal_current_,
            portal_prev_};
  }
  const std::set<std::pair<int, int>>& GetNoRendering() const {
    return no_rendering_buffer_;
  }

 private:
  RectVect domain;
  RectVect resize_queue_;
  blocks Blocks;
  Scal t;
  Scal dt;
  Vect gravity_;
  void calc_forces(size_t i);
  void RHS_bonds();
  void ApplyFrozen();
  std::vector<std::unique_ptr<env_object>> ENVOBJ;
  std::vector<std::vector<size_t>> block_envobj_;
  Vect force_center;
  bool force_enabled;
  bool force_attractive_ = false;
  bool gravity_enable_ = true;
  std::vector<particle> particle_buffer_;
  blocks blocks_buffer_;
  int bonds_prev_particle_id_;
  bool bonds_enabled_ = false;
  int pick_particle_id_;
  bool pick_enabled_ = false;
  Vect pick_pointer_;
  std::set<std::pair<int, int>> bonds_;
  std::set<int> frozen_; // particle id
  bool freeze_enabled_ = false;
  int freeze_last_id_;
  bool renderer_ready_for_next_ = true;
  bool portal_enabled_ = false;
  bool remove_last_portal_;
  int portal_stage_ = 0;
  bool portal_mouse_moving_ = false;
  Vect portal_begin_;
  Vect portal_current_;
  std::pair<Vect, Vect> portal_prev_;
  std::vector<std::array<Portal, 2>> portals_;
  std::vector<int> particle_to_move_; // 1: to move, 0: otherwise
  void UpdatePortalBlocks(Portal& portal);
  // Particles in the blocks of one portal, with their coordinates relative to
  // the portal line: `lambda` spans the portal from 0 to 1 and `offset` is the
  // signed distance to it. Filled by UpdatePortalCache().
  struct PortalCache {
    std::vector<size_t> block, index; // Location in blocks::BlockData.
    std::vector<Vect> position;
    std::vector<Scal> lambda, offset;
    // One entry per non-empty block of the portal, to skip whole blocks at
    // once: the range [begin,end) of their particles in the arrays above, the
    // bounding box of their positions and the range of their `lambda`.
    struct Group {
      size_t begin, end;
      Vect low, high;
      Scal lambda_min, lambda_max;
    };
    std::vector<Group> groups;
  };
  void UpdatePortalCache(const Portal& portal, PortalCache& cache);
  std::array<PortalCache, 2> portal_cache_;
  std::set<std::pair<int, int>> no_rendering_;
  std::set<std::pair<int, int>> no_rendering_buffer_;
  size_t grid_size_ = 0;
};
