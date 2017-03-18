/*
  GEOMETRY
*/
#pragma once

#include <iostream>
#include <cmath>

const double PI=atan(1.0)*4.0;

class vect
{
public:
  double x,y;
  vect() {;}
  vect(double _x, double _y)
  {
    x=_x; y=_y;
  }
  vect& operator*=(double a)
  {
    x*=a; y*=a;
    return *this;
  }
  vect operator*(double a) const
  {
    vect res=*this;
    res*=a;
    return res;
  }
  vect& operator+=(vect v)
  {
    x+=v.x; y+=v.y;
    return *this;
  }
  vect operator+(vect v) const
  {
    vect res=*this;
    res+=v;
    return res;
  }
  vect& operator-=(vect v)
  {
    x-=v.x; y-=v.y;
    return *this;
  }
  vect operator-(vect v) const
  {
    vect res=*this;
    res-=v;
    return res;
  }
  double dot(vect v) const
  {
    return v.x*x+v.y*y;
  }
  double length() const
  {
    return sqrt(x*x+y*y);
  }
  double dist(vect v) const
  {
    return (v-*this).length();
  }
};

class mindex
{
public:
  int i,j;
  mindex() {;}
  mindex(int _i, int _j) : i(_i), j(_j) {;}
  mindex& operator+=(mindex m)
  {
    i+=m.i; j+=m.j;
    return *this;
  }
  mindex operator+(mindex m) const
  {
    mindex res=*this;
    res+=m;
    return res;
  }
};

// A<=B
class rect_vect
{
public:
  vect A,B;
  rect_vect(vect _A, vect _B) : A(_A), B(_B) {;}
  vect size() const
  {
    return B-A;
  }
};

// A<=B
class rect_mindex
{
public:
  mindex A,B;
  rect_mindex(mindex _A, mindex _B) : A(_A), B(_B) {;}
};

class projection
{
public:
  rect_vect Rmath;
  rect_mindex Rscreen;
  projection(rect_vect _Rmath, rect_mindex _Rscreen) : Rmath(_Rmath), Rscreen(_Rscreen) {;}
  mindex convert(vect p)
  {
    mindex res;
    res.i=Rscreen.A.i+int((p.x-Rmath.A.x)/(Rmath.B.x-Rmath.A.x)*(Rscreen.B.i-Rscreen.A.i));
    res.j=Rscreen.A.j+int((p.y-Rmath.A.y)/(Rmath.B.y-Rmath.A.y)*(Rscreen.B.j-Rscreen.A.j));
    return res;
  }
};

std::ostream& operator<<(std::ostream& out, vect v);

struct rgb
{
  float r,g,b;
  rgb(float _r, float _g, float _b) : r(_r), g(_g), b(_b) {;}
  rgb() {;}
};
