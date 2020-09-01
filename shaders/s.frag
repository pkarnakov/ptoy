#version 140
flat in uvec2 upos;
flat in float colorf;
uniform uint pointSize;
void main() {
  vec2 p = gl_FragCoord.xy - upos;
  float d = dot(p, p) / (pointSize * pointSize);
  gl_FragColor = vec4(colorf, 0, 0, 1 - pow(d, 4));
}
