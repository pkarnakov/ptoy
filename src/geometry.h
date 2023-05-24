#pragma once

#include <array>
#include <cmath>
#include <iosfwd>
#include <limits>

using Scal = float;

const Scal PI = atan(1.0) * 4.0;

class Vect {
 public:
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
  friend Vect operator*(Scal k, const Vect& v) {
    return v * k;
  }
  friend std::ostream& operator<<(std::ostream&, Vect);

 public:
  Scal x, y;
};

class MIdx {
 public:
  int i, j;
  MIdx() {
    ;
  }
  MIdx(int _i, int _j) : i(_i), j(_j) {
    ;
  }
  MIdx& operator+=(MIdx m) {
    i += m.i;
    j += m.j;
    return *this;
  }
  MIdx operator+(MIdx m) const {
    MIdx res = *this;
    res += m;
    return res;
  }
};

// A<=B
class RectVect {
 public:
  Vect A, B;
  RectVect(Vect _A, Vect _B) : A(_A), B(_B) {
    ;
  }
  RectVect() = default;
  Vect size() const {
    return B - A;
  }
  bool operator==(const RectVect& other) const {
    return A == other.A && B == other.B;
  }
  bool operator!=(const RectVect& other) const {
    return !(other == *this);
  }
};

// A<=B
class RectMIdx {
 public:
  MIdx A, B;
  RectMIdx(MIdx _A, MIdx _B) : A(_A), B(_B) {
    ;
  }
};

struct rgb {
  float r, g, b;
  rgb(float r_, float g_, float b_) : r(r_), g(g_), b(b_) {}
  rgb() = default;
};

struct GetNanHelper {
  static auto Get(Scal*) {
    return std::numeric_limits<Scal>::quiet_NaN();
  }
  static auto Get(Vect*) {
    return Vect(Get((Scal*)nullptr));
  }
};

template <class T>
T GetNan() {
  return GetNanHelper::Get((T*)nullptr);
}
