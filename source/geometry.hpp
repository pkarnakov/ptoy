/*
  GEOMETRY
*/
#pragma once

#include <iostream>
#include <cmath>

//using Scal = float;
using Scal = double;

const Scal PI=atan(1.0)*4.0;

class vect {
public:
  static const vect kNan;
  Scal x,y;
  vect() {;}
  vect(Scal _x, Scal _y)
  {
    x=_x; y=_y;
  }
  vect& operator*=(Scal a)
  {
    x*=a; y*=a;
    return *this;
  }
  vect operator*(Scal a) const
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
  vect& operator*=(vect v)
  {
    x*=v.x; y*=v.y;
    return *this;
  }
  vect operator*(vect v)
  {
    vect res=*this;
    res*=v;
    return res;
  }
  Scal dot(vect v) const
  {
    return v.x*x+v.y*y;
  }
  Scal length() const
  {
    return sqrt(x*x+y*y);
  }
  Scal dist(vect v) const
  {
    return (v-*this).length();
  }
  bool operator==(vect other) const {
    return x == other.x && y == other.y;
  }
  vect GetNormalized() const {
    const Scal len = length();
    return vect(x / len, y / len);
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
  rect_vect() = default;
  vect size() const
  {
    return B-A;
  }
  bool operator==(const rect_vect& other) const {
    return A == other.A && B == other.B;
  }
  bool operator!=(const rect_vect& other) const {
    return !(other == *this);
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
