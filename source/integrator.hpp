/*
  INTEGRATOR
*/
#pragma once

#include "geometry.hpp"
#include <vector>
#include <iostream>
#include <functional>

using std::endl;
using std::vector;
using std::size_t;

typedef std::function<void(const vector<vect>& X, const vector<vect>& V, double t, vector<vect>& F)> func_RHS;

class integrator
{
protected:
  vector<vect> X, V;
  vector<vect> F;
  func_RHS RHS;
  double t;
  double dt;
public:
  integrator()
  {

  }
  virtual ~integrator()
  {

  }
  virtual void step()
  {
    t+=dt;
  }
  virtual void set_data(const vector<vect>& _X, const vector<vect>& _V)
  {
    X=_X;
    V=_V;
  }
  virtual void get_data(vector<vect>& _X, vector<vect>& _V)
  {
    _X=X;
    _V=V;
  }
  virtual void set_t(double _t)
  {
    t=_t;
  }
  virtual void set_dt(double _dt)
  {
    dt=_dt;
  }
  virtual void set_RHS(func_RHS _RHS)
  {
    RHS=_RHS;
  }
  double get_dt()
  {
    return dt;
  }
  double get_t()
  {
    return t;
  }
  void status(std::ostream& out)
  {
    out<<"Integrator status t="<<t<<", dt="<<dt<<", N="<<X.size()<<endl;
    for(size_t i=0; i<X.size(); ++i)
    {
      out<<i<<":p="<<X[i]<<":v="<<V[i]<<endl;
    }
    out<<endl;
  }
};

class integrator_Euler : public integrator
{
public:
  integrator_Euler()
  {

  }
  ~integrator_Euler()
  {

  }
  void step()
  {
    integrator::step();

    for(size_t i=0; i<X.size(); ++i)
    {
      X[i]+=V[i]*dt;
    }

    RHS(X,V,t,F);

    for(size_t i=0; i<X.size(); ++i)
    {
      V[i]+=(F[i]+V[i]*(-0.1))*dt;
      double limit=100.0;
      if(V[i].length()>limit)
      {
        V[i]=V[i]*(limit/V[i].length());
      }
    }
  }
  void set_data(const vector<vect>& X, const vector<vect>& V)
  {
    integrator::set_data(X, V);
  }
  void get_data(vector<vect>& X, vector<vect>& V)
  {
    integrator::get_data(X, V);
  }
  void set_t(double _t)
  {
    integrator::set_t(_t);
  }
  void set_dt(double _dt)
  {
    integrator::set_dt(_dt);
  }
  void set_RHS(func_RHS _RHS)
  {
    integrator::set_RHS(_RHS);
  }
};
