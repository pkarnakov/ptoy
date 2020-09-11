#version 330 core

uniform sampler2D tex;

out vec4 fragcolor;

void main() {
  fragcolor = texture(tex, gl_FragCoord.xy / 30.);
  fragcolor.a = 0.25;
}
