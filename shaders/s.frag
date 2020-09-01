#version 140
flat in uvec2 upos;
uniform uint pointSize;
void main() {
  vec2 p = gl_FragCoord.xy - upos;
  float d = dot(p, p) / (pointSize * pointSize);
  gl_FragColor = vec4(1, 0, 0, 1 - pow(d, 8));
}
