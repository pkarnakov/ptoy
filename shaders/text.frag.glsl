#version 330 core

uniform uvec2 screenSize;
uniform uvec2 font_size;
uniform float font_offset;
uniform sampler2D tex;

in GEOM {
  flat vec2 lowcorner;
  flat float size;
  flat vec4 color;
  flat int char;
} geom;

out vec4 fragcolor;

void main() {
  vec2 p = (gl_FragCoord.xy - geom.lowcorner) / geom.size;
  if (p.x / font_offset > 0 && p.x / font_offset < 1 && p.y > 0 &&
          p.y < font_size.y || true) {
    float tx = (p.x + geom.char * font_offset) / font_size.x;
    float ty = 1 - p.y / font_size.y;
    fragcolor = geom.color;
    fragcolor.a = 1 - texture(tex, vec2(tx, ty)).x;
  }
}
