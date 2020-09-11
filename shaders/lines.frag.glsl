#version 330 core

in GEOM {
  flat float color;
} geom;

out vec4 fragcolor;

void main() {
  fragcolor = vec4(geom.color, 0, 0, 1);
}
