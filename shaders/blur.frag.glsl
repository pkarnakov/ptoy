#version 330 core

uniform sampler2D tex;
uniform uvec2 screenSize;

out vec4 fragcolor;

void main() {
  vec4 c = texture(tex, (gl_FragCoord.xy) / screenSize);
  float d = 1.5;
  vec4 cxm = texture(tex, (gl_FragCoord.xy + vec2(-d, 0)) / screenSize);
  vec4 cxp = texture(tex, (gl_FragCoord.xy + vec2(d, 0)) / screenSize);
  vec4 cym = texture(tex, (gl_FragCoord.xy + vec2(0, -d)) / screenSize);
  vec4 cyp = texture(tex, (gl_FragCoord.xy + vec2(0, d)) / screenSize);

  float q = 2;
  vec4 cb = (q * c + cxm + cxp + cym + cyp) / (q + 4);
  fragcolor = cb * 0.8;
  fragcolor.a = 0.8;
}
