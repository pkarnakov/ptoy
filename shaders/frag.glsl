#version 330 core

uniform uint pointSize;
uniform sampler2D myTex;

flat in uvec2 upos;
flat in float colorf;

out vec4 fragColor;

void main() {
  vec2 p = gl_FragCoord.xy - upos;
  float d = dot(p, p) / (pointSize * pointSize);
  vec4 c = vec4(colorf, 0, 0, 1 - pow(d, 4));
  c.y = texture(myTex, (p / pointSize + 1) * 0.5).x * colorf;
  fragColor = c;
}
