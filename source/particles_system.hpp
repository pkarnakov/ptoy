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

#include <x86intrin.h>

using std::endl;
using std::vector;
using std::size_t;
using std::min;

const Scal kRadius = 0.02;
const Scal kSigma = 0.5;
const Scal kMass = 0.01;

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
  Scal eps;
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
  line(vect _A, vect _B, Scal _eps) : A(_A), B(_B), eps(_eps) {;}
  vect F(vect p, vect /*v*/, Scal R, Scal sigma) override {
    return F12(p, vect(0.,0.), GetNearest(p), vect(0.,0.), sigma, R);
  }
  bool IsClose(vect p, Scal R) {
    return GetNearest(p).dist(p) < R + kRadius;
  }

};

class particles_system
{
 public:
  particles_system(); 
  ~particles_system();
  std::vector<particle> GetParticles();
  void AddEnvObj(env_object* env);
  void ClearEnvObj() { ENVOBJ.clear(); }
  void UpdateEnvObj();
  void SetDomain(rect_vect new_domain) { 
    std::lock_guard<std::mutex> lg(m_step);
    domain = new_domain; 
    // TODO: window resize
    //Blocks.SetDomain(domain);
  }
  void status(std::ostream& out);
  void step(Scal time_target);
  void SetForce(vect center, bool enabled);
  void SetForce(vect center);
  void SetForce(bool enabled);
  Scal GetTime() const { return t; }
  rect_vect GetDomain() const { return domain; }
  mutable std::mutex m_ENVOBJ;
  mutable std::mutex m_step;
  size_t GetNumParticles() const { 
    //std::lock_guard<std::mutex> lg(m_step);
    return Blocks.GetNumParticles(); 
  }
  size_t GetNumPerCell() const {
    return Blocks.GetNumPerCell();
  }

 private:
  rect_vect domain;
  blocks Blocks;
  Scal t;
  Scal dt;
  vect g;
  void RHS(size_t i);
  vector<std::unique_ptr<env_object>> ENVOBJ;
  std::vector<std::vector<size_t>> block_envobj_;
  vect force_center;
  bool force_enabled;
};

