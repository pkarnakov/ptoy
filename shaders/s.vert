#version 140
in vec2 point;
uniform uvec2 screenSize;
flat out uvec2 upos;
uniform uint pointSize;
void main() {
  float x = point.x;
  float y = point.y;
  upos.x = uint((x + 1) * 0.5 * screenSize.x);
  upos.y = uint((y + 1) * 0.5 * screenSize.y);
  gl_Position = vec4(x, y, 0, 1);
  gl_PointSize = float(pointSize) * 2;
}
