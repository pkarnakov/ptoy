#version 330 core

uniform sampler2D tex;

out vec4 fragcolor;

void main() {
  fragcolor = texture(tex, gl_FragCoord.xy / 800 - vec2(0, 0));
  fragcolor.a = 0.6;
}
