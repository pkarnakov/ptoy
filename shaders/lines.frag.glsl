#version 330 core

in float geomcolor;

out vec4 fragcolor;

void main() {
  fragcolor = vec4(geomcolor, 0, 0, 1);
}
