#version 330

uniform uvec2 screenSize;
uniform uvec2 font_size;
uniform float font_offset;

layout (points) in;
in VERT {
  flat vec2 lowcorner;
  flat float size;
  flat vec4 color;
  flat int char;
} vert[];

layout  (triangle_strip, max_vertices = 4) out;

out GEOM {
  flat vec2 lowcorner;
  flat float size;
  flat vec4 color;
  flat int char;
} geom;

void main() {
  vec4 p = gl_in[0].gl_Position;
  vec2 s = vec2(font_offset, font_size.y) / screenSize * 2 * vert[0].size;

  geom.lowcorner = vert[0].lowcorner;
  geom.size = vert[0].size;
  geom.color = vert[0].color;
  geom.char = vert[0].char;

  gl_Position = p + vec4(0, 0, 0, 0);
  EmitVertex();

  gl_Position = p + vec4(s.x, 0, 0, 0);
  EmitVertex();

  gl_Position = p + vec4(0, s.y, 0, 0);
  EmitVertex();

  gl_Position = p + vec4(s.x, s.y, 0, 0);
  EmitVertex();

  EndPrimitive();
}
