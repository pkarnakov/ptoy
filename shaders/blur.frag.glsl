#version 330 core

uniform sampler2D tex;

out vec4 fragcolor;

void main() {
  vec4 c = texture(tex, gl_FragCoord.xy / 800 - vec2(0, 0));
  fragcolor = c * 0.5;
  fragcolor.a = 0.9;
  fragcolor.x += c.x * 0.2;
  fragcolor.y += c.x * 0.4;
  fragcolor.z += c.x * 0.7;
}
