#pragma once

#include <array>
#include <cmath>
#include <iostream>

using Scal = float;

const Scal PI = atan(1.0) * 4.0;

class Vect {
 public:
  static const Vect kNan;
  Scal x, y;
  Vect() {}
  Vect(Scal x_) : x(x_), y(x_) {}
  Vect(Scal x_, Scal y_) : x(x_), y(y_) {}
  Vect& operator*=(Scal a) {
    x *= a;
    y *= a;
    return *this;
  }
  Vect operator*(Scal a) const {
    Vect res(*this);
    res *= a;
    return res;
  }
  Vect& operator+=(Vect v) {
    x += v.x;
    y += v.y;
    return *this;
  }
  Vect operator+(Vect v) const {
    Vect res(*this);
    res += v;
    return res;
  }
  Vect& operator-=(Vect v) {
    x -= v.x;
    y -= v.y;
    return *this;
  }
  Vect operator-(Vect v) const {
    Vect res(*this);
    res -= v;
    return res;
  }
  Vect& operator*=(Vect v) {
    x *= v.x;
    y *= v.y;
    return *this;
  }
  Vect operator*(Vect v) const {
    Vect res(*this);
    res *= v;
    return res;
  }
  Scal dot(Vect v) const {
    return v.x * x + v.y * y;
  }
  Scal length() const {
    return sqrt(x * x + y * y);
  }
  Scal dist(Vect v) const {
    return (v - *this).length();
  }
  bool operator==(Vect other) const {
    return x == other.x && y == other.y;
  }
  Vect GetNormalized() const {
    const Scal len = length();
    return Vect(x / len, y / len);
  }
  template <class T>
  operator std::array<T, 2>() {
    return {x, y};
  }
};

class mindex {
 public:
  int i, j;
  mindex() {
    ;
  }
  mindex(int _i, int _j) : i(_i), j(_j) {
    ;
  }
  mindex& operator+=(mindex m) {
    i += m.i;
    j += m.j;
    return *this;
  }
  mindex operator+(mindex m) const {
    mindex res = *this;
    res += m;
    return res;
  }
};

// A<=B
class rect_vect {
 public:
  Vect A, B;
  rect_vect(Vect _A, Vect _B) : A(_A), B(_B) {
    ;
  }
  rect_vect() = default;
  Vect size() const {
    return B - A;
  }
  bool operator==(const rect_vect& other) const {
    return A == other.A && B == other.B;
  }
  bool operator!=(const rect_vect& other) const {
    return !(other == *this);
  }
};

// A<=B
class rect_mindex {
 public:
  mindex A, B;
  rect_mindex(mindex _A, mindex _B) : A(_A), B(_B) {
    ;
  }
};

std::ostream& operator<<(std::ostream& out, Vect v);

struct rgb {
  float r, g, b;
  rgb(float _r, float _g, float _b) : r(_r), g(_g), b(_b) {
    ;
  }
  rgb() {
    ;
  }
};
