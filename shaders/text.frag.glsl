#version 330 core

uniform uvec2 screenSize;
uniform uvec2 font_size;
uniform sampler2D tex;

in GEOM {
  flat vec2 lowcorner;
  flat float size;
  flat vec4 color;
  flat float char;
  flat float width;
} geom;

out vec4 fragcolor;

void main() {
  vec2 p = (gl_FragCoord.xy - geom.lowcorner) / geom.size;
  if (p.x / geom.width > 0 && p.x / geom.width < 1 && p.y > 0 &&
          p.y < font_size.y || true) {
    float tx = (p.x + geom.char) / font_size.x;
    float ty = 1 - p.y / font_size.y;
    fragcolor = geom.color;
    fragcolor.a = 1 - texture(tex, vec2(tx, ty)).x;
  }
}
