#include <iostream>

#include "geometry.hpp"

std::ostream& operator<<(std::ostream& out, Vect v) {
  out << "(" << v.x << "," << v.y << ")";
  return out;
}
