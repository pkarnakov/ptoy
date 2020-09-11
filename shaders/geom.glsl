#version 330

layout (points) in;
flat in uvec2 upos[];
flat in float colorf[];

layout  (points, max_vertices = 1) out;
flat out uvec2 geom_upos;
flat out float geom_colorf;

void main() {
  gl_Position = gl_in[0].gl_Position;
  EmitVertex();
  EndPrimitive();
}
