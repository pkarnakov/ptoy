#version 140
in vec2 LVertexPos2D;
void main() {
  float x = LVertexPos2D.x;
  float y = LVertexPos2D.y;
  gl_Position = vec4(x, x > 0 ? y : y - 0.3, 0, 1);
}
