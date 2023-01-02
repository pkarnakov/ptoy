#include "geometry.hpp"

std::ostream& operator<<(std::ostream& out, Vect v) {
  out << "(" << v.x << "," << v.y << ")";
  return out;
}

const Scal ScalNan = std::nan("");
const Vect Vect::kNan = Vect(ScalNan, ScalNan);
