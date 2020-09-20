#version 330 core

in GEOM {
  flat vec4 color;
} geom;

out vec4 fragcolor;

void main() {
  fragcolor = geom.color;
}
