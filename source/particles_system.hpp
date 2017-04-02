/*
  PARTICLES SYSTEM
*/
#pragma once

#include "geometry.hpp"
#include <vector>
#include <iostream>
#include <functional>
#include <memory>
#include <mutex>
#include "blocks.hpp"
#include <atomic>
#include <set>

#include <x86intrin.h>

using std::endl;
using std::vector;
using std::size_t;
using std::min;

const Scal kRadius = 0.02;
const Scal kSigma = 0.5;
const Scal kMass = kRadius * kRadius * 100.;
const Scal kPointForce = 0.1;
const Scal kPointForceAttractive = 1.1;
const Scal kDissipation = 0.01;
const Scal kTimeStep = 0.0003;
const Scal kBlockSize = 4. * kRadius;
const Scal kGravity = 10.;

const int kParticleIdNone = -1;

template<class T>
T sqr(T a)
{
  return a*a;
}

void TestUni();

struct alignas(64) particle
{
  vect p;
  vect v;
  vect f;
  vect p0;
  vect v0;
  //unsigned int layers_mask;
  //rgb color;
  particle() {;}
  particle(vect _p, vect _v, Scal /*_m*/, Scal /*_r*/, 
      Scal /*_sigma*/, 
      unsigned int /*_layers_mask*/, 
      rgb /*_color*/)
    : p(_p), v(_v)//, layers_mask(_layers_mask), color(_color)
  {;}
};

vect F12(vect p1, vect v1, vect p2, vect v2, Scal sigma, Scal R);

class env_object
{
public:
  virtual vect F(vect p, vect v, Scal R, Scal sigma) = 0;
  virtual bool IsClose(vect p, Scal R) = 0;
};

class line : public env_object
{
  vect A, B;
  //Scal eps;
  vect GetNearest(vect p) {
    vect Q;
    Scal lambda = (B - A).dot(p - A) / (B - A).dot(B - A);
    if (lambda > 0. && lambda < 1.) {
      Q = A + (B-A) * lambda;
    } else {
      Q = (p.dist(A) < p.dist(B)) ? A : B;
    }
    return Q;
  }
public:
  line(vect _A, vect _B) : A(_A), B(_B) {;}
  vect F(vect p, vect /*v*/, Scal R, Scal sigma) override {
    return F12(p, vect(0.,0.), GetNearest(p), vect(0.,0.), sigma, R);
  }
  bool IsClose(vect p, Scal R) override {
    return GetNearest(p).dist(p) < R + kRadius;
  }

};

class particles_system
{
 public:
  particles_system(); 
  ~particles_system();
  void SetParticleBuffer();
  void AddEnvObj(env_object* env);
  void ClearEnvObj() { ENVOBJ.clear(); }
  void UpdateEnvObj();
  void ResetEnvObjFrame(rect_vect new_domain) {
    const vect A = new_domain.A, B = new_domain.B;
    ClearEnvObj();
    //std::cout << "envobj:" << domain.A << " " << domain.B << std::endl;
    AddEnvObj(new line(vect(A.x, A.y), vect(B.x, A.y)));
    AddEnvObj(new line(vect(A.x, B.y), vect(B.x, B.y)));
    AddEnvObj(new line(vect(A.x, A.y), vect(A.x, B.y)));
    AddEnvObj(new line(vect(B.x, A.y), vect(B.x, B.y)));
    UpdateEnvObj();
  }
  void SetDomain(rect_vect new_domain) { 
    domain = new_domain;
    Blocks.SetDomain(domain);
  }
  void PushResize(rect_vect new_domain) {
    resize_queue_ = new_domain;
  }
  const std::vector<std::pair<int, int>>& GetBonds() const {
    return bonds_;
  }
  const std::set<int>& GetFrozen() const {
    return frozen_;
  }
  void status(std::ostream& out);
  void step(Scal time_target, const std::atomic<bool>& quit);
  void SetForce(vect center, bool enabled);
  void SetForce(vect center);
  void SetForce(bool enabled);
  void SetForceAttractive(bool value) {
    force_attractive_ = value;
  }
  void BondsStart(vect point);
  void BondsMove(vect point);
  void BondsStop(vect point);
  void FreezeStart(vect point);
  void FreezeMove(vect point);
  void FreezeStop(vect point);
  void PickStart(vect point);
  void PickMove(vect point);
  void PickStop(vect point);
  Scal GetTime() const { return t; }
  size_t GetNumSteps() const { return static_cast<size_t>(t / dt); } 
  rect_vect GetDomain() const { return domain; }

  // These are only for external use (TODO: check or ensure)
  bool GetGravity() const {
    return gravity_enable_;
  }
  void InvertGravity() {
    gravity_enable_ = !gravity_enable_;
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
  const std::vector<std::pair<size_t, size_t>> GetBlockById() const {
    return blocks_buffer_.GetBlockById();
  }
  const blocks::BlockData& GetBlockData() const {
    return blocks_buffer_.GetData();
  };
  void SetRendererReadyForNext(bool value) {
    renderer_ready_for_next_ = value;
  }

  mutable std::mutex m_buffer_;
 private:
  rect_vect domain;
  rect_vect resize_queue_;
  blocks Blocks;
  Scal t;
  Scal dt;
  vect g;
  void RHS(size_t i);
  void RHS_bonds();
  vector<std::unique_ptr<env_object>> ENVOBJ;
  std::vector<std::vector<size_t>> block_envobj_;
  vect force_center;
  bool force_enabled;
  bool force_attractive_ = false;
  bool gravity_enable_ = true;
  std::vector<particle> particle_buffer_;
  blocks blocks_buffer_;
  int bonds_prev_particle_id_;
  bool bonds_enabled_ = false;
  int pick_particle_id_;
  bool pick_enabled_ = false;
  vect pick_pointer_;
  std::vector<std::pair<int, int>> bonds_;
  std::set<int> frozen_; // particle id
  bool freeze_enabled_ = false;
  int freeze_last_id_;
  bool renderer_ready_for_next_ = true;
};

