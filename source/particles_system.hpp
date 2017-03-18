/*
  PARTICLES SYSTEM
*/
#pragma once

#include "geometry.hpp"
#include "integrator.hpp"
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
  vect F(vect p, vect v, double R, double sigma)
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
  vector<vect> X, V;
  blocks Blocks;
  double t;
  double dt;
  vect g;
  void transfer_data(const vector<particle>& P, vector<vect>& X, vector<vect>& V)
  {
    size_t N=P.size();
    X.resize(N);
    V.resize(N);
    for(size_t i=0; i<N; ++i)
    {
      X[i]=P[i].p;
      V[i]=P[i].v;
    }
  }

public:
  std::unique_ptr<integrator> INT;
  vector<particle> P;
  vector<std::unique_ptr<env_object>> ENVOBJ;
  particles_system() : Blocks(rect_vect(vect(-1.,-1.),vect(1.,1.)), vect(0.2, 0.2), 100)
  {
    INT=std::unique_ptr<integrator>(new integrator_Euler);

    // place particles in the domain
    double r=0.02;
    int N=100;
    for(int i=0; i<N; ++i)
    {
      int row=11;
      P.push_back(particle(vect((i%row*2.0+1.0)*r-1.0, (i/row*2.0+1.0)*r-1.0), vect(0.1, 0.0), 0.01, r, 10000.0, 1+i%2*0, rgb(1.0, i%2, 0.0)));
    }

    transfer_data(P, X,V);

    t=0.0;
    dt=0.0001;
    g=vect(0.0, -10.0);

    INT->set_data(X, V);
    INT->set_t(t);
    INT->set_dt(dt);

    using std::placeholders::_1;
    using std::placeholders::_2;
    using std::placeholders::_3;
    using std::placeholders::_4;
    INT->set_RHS(std::bind(&particles_system::RHS, this, _1, _2, _3, _4));
  }
  ~particles_system()
  {

  }
  void status(std::ostream& out)
  {
    out<<"Particles system"<<std::endl<<"Particles number = "<<P.size()<<std::endl;
  }
  void step()
  {
    INT->step();

    INT->get_data(X,V);

    for(std::size_t i=0; i<P.size(); ++i)
    {
      P[i].p=X[i];
      P[i].v=V[i];
    }

    Blocks.arrange(X);

    t=INT->get_t();
  }
  void RHS(const vector<vect>& X, const vector<vect>& V, double t, vector<vect>& F) const
  {
    F.resize(X.size());

    for(size_t i=0; i<F.size(); ++i)
    {
      F[i]=g*P[i].m;
      for(size_t k=0; k<ENVOBJ.size(); ++k)
      {
        auto& obj=ENVOBJ[k];
        F[i]+=obj->F(X[i],V[i],P[i].r,P[i].sigma);
      }
    }

    for(int bj=0; bj<Blocks.B.msize().j; ++bj)
    for(int bi=0; bi<Blocks.B.msize().i; ++bi)
    {
      mindex m(bi,bj);
      const std::vector<int>& b1=Blocks.B[m];
      for(int w1=0; w1<b1.size(); ++w1)
      {
        int p1=b1[w1]; // first particle

        for(size_t k=0; k<Blocks.NEAR.size(); ++k)
        {
          mindex mnear=m+Blocks.NEAR[k];
          if(Blocks.B.valid(mnear))
          {
            const std::vector<int>& b2=Blocks.B[mnear];
            for(int w2=0; w2<b2.size(); ++w2)
            {
              int p2=b2[w2]; // second particle
              if(p1!=p2 && (P[p1].layers_mask & P[p2].layers_mask))
              {
                F[p1]+=F12(X[p1], V[p1], X[p2], V[p2], 0.5*(P[p1].sigma+P[p2].sigma), P[p1].r+P[p2].r);
              }
            }
          }
        }
      }
    }

    for(size_t i=0; i<F.size(); ++i)
    {
      F[i]*=1.0/P[i].m;
    }
  }
};
