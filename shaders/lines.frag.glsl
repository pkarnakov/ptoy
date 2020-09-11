#version 330 core

flat in float fragcolor;

out vec4 FragColor;

void main() {
  FragColor = vec4(fragcolor, 0, 0, 1);
}
