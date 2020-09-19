#version 330 core

uniform sampler2D tex;
uniform uvec2 screenSize;

out vec4 fragcolor;

void main() {
  vec4 c = texture(tex, gl_FragCoord.xy / screenSize);
  fragcolor = c * 0.5;
  fragcolor.a = 0.9;
  fragcolor.x += c.x * 0.2;
  fragcolor.y += c.x * 0.4;
  fragcolor.z += c.x * 0.7;
}
