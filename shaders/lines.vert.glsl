#version 330

uniform vec2 domain0;
uniform vec2 domain1;

in vec2 point;
in vec4 color;
in float width;

out VERT {
  flat vec4 color;
  flat float width;
} vert;

void main() {
  vec2 rel = (point - domain0) / (domain1 - domain0);
  gl_Position = vec4(rel * 2 - 1, 0, 1);
  vert.color = color;
  vert.width = width;
}
