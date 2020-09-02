#version 140

uniform uvec2 screenSize;
uniform uint pointSize;
uniform vec2 domain0;
uniform vec2 domain1;

in vec2 point;
in float color;

flat out uvec2 upos;
flat out float colorf;
void main() {
  vec2 rel = (point - domain0) / (domain1 - domain0);
  upos = uvec2(rel * screenSize);
  gl_Position = vec4(rel * 2 - 1, 0, 1);
  gl_PointSize = float(pointSize) * 2 + 1;
  colorf = color;
}
