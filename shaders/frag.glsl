#version 330 core

uniform uint pointSize;
uniform sampler2D tex;
uniform vec3 facecolor;

flat in uvec2 geom_upos;
flat in float geom_colorf;

out vec4 fragColor;

void main() {
  vec2 p = gl_FragCoord.xy - geom_upos;
  float d = dot(p, p) / (pointSize * pointSize);
  float q = geom_colorf;
  vec4 c =
      vec4(q * facecolor.x, q * facecolor.y, q * facecolor.z, 1 - pow(d, 2));
  float t = texture(tex, (p / pointSize + 1) * 0.5).x;
  c *= 0.7 + t;
  fragColor = c;
}
