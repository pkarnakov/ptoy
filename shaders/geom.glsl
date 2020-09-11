#version 330

// FIXME: this is an identity shader, it copies input to output.
// Without it, 'lines.geom.glsl' is not executed while drawing lines.

layout (points) in;
flat in uvec2 upos[];
flat in float colorf[];

layout  (points, max_vertices = 1) out;
flat out uvec2 geom_upos;
flat out float geom_colorf;

void main() {
  gl_Position = gl_in[0].gl_Position;
  gl_PointSize = gl_in[0].gl_PointSize;
  geom_upos = upos[0];
  geom_colorf = colorf[0];
  EmitVertex();
  EndPrimitive();
}
