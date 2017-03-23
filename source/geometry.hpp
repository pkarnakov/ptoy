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


std::ostream& operator<<(std::ostream& out, vect v);

struct rgb
{
  float r,g,b;
  rgb(float _r, float _g, float _b) : r(_r), g(_g), b(_b) {;}
  rgb() {;}
};
