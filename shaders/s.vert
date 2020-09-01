#version 140

uniform uvec2 screenSize;
uniform uint pointSize;

in vec2 point;
in float color;

flat out uvec2 upos;
flat out float colorf;
void main() {
  float x = point.x;
  float y = point.y;
  upos.x = uint((x + 1) * 0.5 * screenSize.x);
  upos.y = uint((y + 1) * 0.5 * screenSize.y);
  gl_Position = vec4(x, y, 0, 1);
  gl_PointSize = float(pointSize) * 2;
  colorf = color;
}
