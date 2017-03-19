/*
  PARTICLES SYSTEM
*/
#pragma once

#include "geometry.hpp"
#include <vector>
#include <iostream>
#include <functional>
#include "blocks.hpp"
#include <memory>


using std::endl;
using std::vector;
using std::size_t;
using std::min;

template<class T>
T sqr(T a)
{
  return a*a;
}

struct particle
{
  vect p;
  vect v;
  double m;
  double r;
  double sigma;
  unsigned int layers_mask;
  rgb color;
  particle() {;}
  particle(vect _p, vect _v, double _m, double _r, double _sigma, unsigned int _layers_mask, rgb _color)
    : p(_p), v(_v), m(_m), r(_r), sigma(_sigma), layers_mask(_layers_mask) , color(_color)
  {;}
};

vect F12(vect p1, vect v1, vect p2, vect v2, double sigma, double R);

class env_object
{
public:
  virtual vect F(vect p, vect v, double R, double sigma) = 0;
};

class line : public env_object
{
  vect A, B;
  double eps;
public:
  line(vect _A, vect _B, double _eps) : A(_A), B(_B), eps(_eps) {;}
  vect F(vect p, vect /*v*/, double R, double sigma)
  {
    vect Q;
    double lambda=(B-A).dot(p-A)/(B-A).dot(B-A);
    if(lambda>0. && lambda<1.)
    {
      Q=A+(B-A)*lambda;
    }
    else
    {
      Q=(p.dist(A)<p.dist(B))?A:B;
    }

    return F12(p, vect(0.,0.), Q, vect(0.,0.), sigma, R);
  }
};

class particles_system
{
 public:
  particles_system(); 
  ~particles_system();
  std::vector<particle> GetParticles() const;
  void AddEnvObj(env_object* env);
  void status(std::ostream& out);
  void step();
  void SetForce(vect center, bool enabled);
  void SetForce(vect center);
  void SetForce(bool enabled);

 private:
  vector<vect> X, V;
  blocks Blocks;
  double t;
  double dt;
  vect g;
  void transfer_data(const vector<particle>& P, vector<vect>& X, vector<vect>& V);
  std::vector<vect> RHS() const;
  vector<std::unique_ptr<env_object>> ENVOBJ;
  vector<particle> P;
  vect force_center;
  bool force_enabled;
};


