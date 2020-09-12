#version 330 core

uniform uint pointSize;
uniform sampler2D myTex;

flat in uvec2 geom_upos;
flat in float geom_colorf;

out vec4 fragColor;

void main() {
  vec2 p = gl_FragCoord.xy - geom_upos;
  float d = dot(p, p) / (pointSize * pointSize);
  float q = geom_colorf;
  vec4 c = vec4(q, q * 0.5, q * 0.0, 1 - pow(d, 2));
  //c *= 0.5 + texture(myTex, (p / pointSize + 1) * 0.5);
  fragColor = c;
}
