#version 330

uniform vec2 domain0;
uniform vec2 domain1;

layout (lines) in;
in float vertcolor[];

layout  (triangle_strip, max_vertices = 4) out;
out float geomcolor;

void main() {
  vec4 d = gl_in[1].gl_Position - gl_in[0].gl_Position;
  float lx = domain1.x - domain0.x;
  float ly = domain1.y - domain0.y;
  vec4 n = normalize(vec4(d.y * ly, -d.x * lx, 0, 0)) * 0.01;
  n.x /= lx;
  n.y /= ly;

  gl_Position = gl_in[0].gl_Position - n;
  geomcolor = vertcolor[0];
  EmitVertex();

  gl_Position = gl_in[1].gl_Position - n;
  geomcolor = vertcolor[1];
  EmitVertex();


  gl_Position = gl_in[0].gl_Position + n;
  geomcolor = vertcolor[0];
  EmitVertex();

  gl_Position = gl_in[1].gl_Position + n;
  geomcolor = vertcolor[1];
  EmitVertex();

  EndPrimitive();
}
