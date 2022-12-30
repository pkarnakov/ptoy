#version 330 core

uniform uvec2 screenSize;
uniform uvec2 font_size;
uniform sampler2D tex;

in GEOM {
  flat vec2 lowcorner;
  flat vec2 size;
  flat vec4 color;
} geom;
out vec4 fragcolor;

void main() {
  vec2 r = geom.size * 0.5;
  vec2 p = gl_FragCoord.xy - geom.lowcorner;
  vec2 d = abs(p - r);
  float q = 4;
  float rs = 6;
  float u = pow(max(0, d.x - r.x + rs), 2) + pow(max(0, d.y - r.y + rs), 2);
  u = 1 - u / (rs * rs);
  // u=0 is boundary of smoothed rectangle, positive inside

  fragcolor = geom.color;
  fragcolor.a = u;
}
