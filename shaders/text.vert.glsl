#version 330

uniform uvec2 screenSize;

in vec2 lowcorner;
in float size;
in vec4 color;
in int char;

out VERT {
  flat vec2 lowcorner;
  flat float size;
  flat vec4 color;
  flat int char;
} vert;

void main() {
  gl_Position = vec4(lowcorner / screenSize * 2 - 1, 0, 1);
  vert.lowcorner = lowcorner;
  vert.size = size;
  vert.color = color;
  vert.char = char;
}
