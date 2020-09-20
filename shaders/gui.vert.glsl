#version 330

uniform uvec2 screenSize;

in vec2 lowcorner;
in vec2 size;
in vec4 color;

out VERT {
  flat vec2 lowcorner;
  flat vec2 size;
  flat vec4 color;
} vert;

void main() {
  gl_Position = vec4(lowcorner / screenSize * 2 - 1, 0, 1);
  vert.lowcorner = lowcorner;
  vert.size = size;
  vert.color = color;
}
