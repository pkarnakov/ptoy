#version 330

uniform vec2 domain0;
uniform vec2 domain1;

in vec2 point;
in float color;
in float width;

flat out float fragcolor;

void main() {
  vec2 rel = (point - domain0) / (domain1 - domain0);
  gl_Position = vec4(rel * 2 - 1, 0, 1);
  fragcolor = width + color;
}
