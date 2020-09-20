#version 330

uniform uvec2 screenSize;

layout (points) in;
in VERT {
  flat vec2 lowcorner;
  flat vec2 size;
  flat vec4 color;
} vert[];

layout  (triangle_strip, max_vertices = 4) out;

out GEOM {
  flat vec2 lowcorner;
  flat vec2 size;
  flat vec4 color;
} geom;

void main() {
  vec4 p = gl_in[0].gl_Position;
  vec2 s = vert[0].size / screenSize * 2;

  geom.lowcorner = vert[0].lowcorner;
  geom.size = vert[0].size;
  geom.color = vert[0].color;

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
