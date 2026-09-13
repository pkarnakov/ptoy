var output = document.getElementById('output');
var outputerr = document.getElementById('outputerr');
var g_tmp_canvas;

// Particles.
var GetParticles;
var g_particles;
var g_particles_ptr;
var g_particles_max_size = 10000;

// Portals.
var GetPortals;
var g_portals;
var g_portals_ptr;
var g_portals_max_size = 10000;

// Bonds.
var GetBonds;
var g_bonds;
var g_bonds_ptr;
var g_bonds_max_size = 10000;

// Frozen particles.
var GetFrozen;
var g_frozen;
var g_frozen_ptr;
var g_frozen_max_size = 10000;


// Colors.
var c_red = "#ff1f5b";
var c_green = "#00cd6c";
var c_blue = "#009ade";
var c_orange = "#f28522";
var c_gray = "#a0b1ba";
var c_white = "#ffffff";
var c_black = "#000000";
var c_background = "#50585d";

// Control.
var SendKeyDown;
var SendMouseMotion;
var SendMouseDown;
var SendMouseUp;

// Gravity.
var GetGravity;
var SetGravity;
var SetGravityVect;
var g_accel;
var flag_accel = false;

// Misc
var SetPause;
var flag_pause = false;
var GetMouseMode;
var Init;

// WebGL renderer. Null until initGl() succeeds, drawing falls back to the 2D
// canvas while it is.
var gl = null;
var g_circle;  // Program for particles, drawn as instanced quads.
var g_line;    // Program for portals and bonds, drawn as instanced quads.
var g_quad;    // Unit quad shared by both, one instance of the geometry.
var g_gl_particles;
var g_gl_portals;
var g_gl_bonds;
var g_gl_frozen;


// "#rrggbb" or "#rrggbbaa" to [r, g, b, a] in 0..1. Memoized, the callers are
// in the draw loop.
var g_colors = {};
function parseColor(hex) {
  if (!g_colors[hex]) {
    const s = hex.replace('#', '');
    const n = parseInt(s, 16);
    g_colors[hex] = s.length == 8
        ? [(n >>> 24) / 255, ((n >>> 16) & 255) / 255,
           ((n >>> 8) & 255) / 255, (n & 255) / 255]
        : [((n >>> 16) & 255) / 255, ((n >>> 8) & 255) / 255,
           (n & 255) / 255, 1];
  }
  return g_colors[hex];
}

// Pixel coordinates are up to 65535, well past what mediump resolves, so both
// stages run in highp.
const kShaderCircleVert = `#version 300 es
in vec2 a_corner;
in vec2 a_center;
uniform vec2 u_resolution;
uniform float u_radius;
out vec2 v_offset;
void main() {
  // One pixel of margin leaves room for the outer half of the antialiased rim.
  v_offset = a_corner * (u_radius + 1.0);
  vec2 c = (a_center + v_offset) / u_resolution * 2.0 - 1.0;
  gl_Position = vec4(c.x, -c.y, 0.0, 1.0);
}
`;

const kShaderCircleFrag = `#version 300 es
precision highp float;
in vec2 v_offset;
uniform float u_radius;
uniform vec4 u_color;
out vec4 o_color;
void main() {
  // Coverage over one pixel centred on the rim, so the solid part of the disc
  // reaches u_radius as the filled arc of the 2D canvas does.
  float a = 1.0 - smoothstep(u_radius - 0.5, u_radius + 0.5, length(v_offset));
  if (a <= 0.0) {
    discard;
  }
  o_color = vec4(u_color.rgb, u_color.a * a);
}
`;

// Lines are expanded to quads here because WebGL caps gl_LineWidth at 1.
const kShaderLineVert = `#version 300 es
in vec2 a_corner;
in vec4 a_segment;
uniform vec2 u_resolution;
uniform float u_width;
flat out int v_instance;
void main() {
  vec2 a = a_segment.xy;
  vec2 d = a_segment.zw - a;
  float len = length(d);
  // An empty segment collapses to a zero-area quad and rasterizes nothing,
  // which is how incomplete portal pairs stay invisible.
  vec2 normal = len > 0.0 ? vec2(-d.y, d.x) / len : vec2(0.0);
  vec2 p = a + d * (a_corner.x + 1.0) * 0.5 + normal * a_corner.y * u_width * 0.5;
  vec2 c = p / u_resolution * 2.0 - 1.0;
  gl_Position = vec4(c.x, -c.y, 0.0, 1.0);
  v_instance = gl_InstanceID;
}
`;

const kShaderLineFrag = `#version 300 es
precision highp float;
flat in int v_instance;
uniform vec4 u_color0;
uniform vec4 u_color1;
out vec4 o_color;
void main() {
  o_color = (v_instance % 2 == 0) ? u_color0 : u_color1;
}
`;

function createProgram(vert, frag, uniforms) {
  let compile = function(type, source) {
    let shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      throw new Error(gl.getShaderInfoLog(shader));
    }
    return shader;
  };
  let prog = gl.createProgram();
  gl.attachShader(prog, compile(gl.VERTEX_SHADER, vert));
  gl.attachShader(prog, compile(gl.FRAGMENT_SHADER, frag));
  gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
    throw new Error(gl.getProgramInfoLog(prog));
  }
  let res = {prog: prog, loc: {}};
  uniforms.forEach(u => { res.loc[u] = gl.getUniformLocation(prog, u); });
  return res;
}

// One instance buffer plus the vertex array binding it to the unit quad.
// `components` is 2 for a point and 4 for a segment.
function createInstances(prog, name, components, max_size) {
  let buffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
  gl.bufferData(gl.ARRAY_BUFFER, max_size * 2, gl.DYNAMIC_DRAW);

  let vao = gl.createVertexArray();
  gl.bindVertexArray(vao);
  let corner = gl.getAttribLocation(prog.prog, 'a_corner');
  gl.bindBuffer(gl.ARRAY_BUFFER, g_quad);
  gl.enableVertexAttribArray(corner);
  gl.vertexAttribPointer(corner, 2, gl.FLOAT, false, 0, 0);
  let instance = gl.getAttribLocation(prog.prog, name);
  gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
  gl.enableVertexAttribArray(instance);
  // The packed pixel coordinates are read as unnormalized floats, so the
  // buffers filled by GetParticles() and friends are used as they are.
  gl.vertexAttribPointer(instance, components, gl.UNSIGNED_SHORT, false, 0, 0);
  gl.vertexAttribDivisor(instance, 1);
  gl.bindVertexArray(null);

  return {buffer: buffer, vao: vao, components: components};
}

function initGl(canvas) {
  gl = canvas.getContext(
      'webgl2', {alpha: false, antialias: true, depth: false, stencil: false});
  if (!gl) {
    return false;
  }
  try {
    g_quad = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, g_quad);
    gl.bufferData(
        gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]),
        gl.STATIC_DRAW);

    g_circle = createProgram(
        kShaderCircleVert, kShaderCircleFrag,
        ['u_resolution', 'u_radius', 'u_color']);
    g_line = createProgram(
        kShaderLineVert, kShaderLineFrag,
        ['u_resolution', 'u_width', 'u_color0', 'u_color1']);

    g_gl_particles =
        createInstances(g_circle, 'a_center', 2, g_particles_max_size);
    g_gl_frozen = createInstances(g_circle, 'a_center', 2, g_frozen_max_size);
    g_gl_portals = createInstances(g_line, 'a_segment', 4, g_portals_max_size);
    g_gl_bonds = createInstances(g_line, 'a_segment', 4, g_bonds_max_size);
  } catch (error) {
    console.log(error);
    gl = null;
    return false;
  }

  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  return true;
}

// Uploads `count` instances from `view` and draws them as quads.
function drawInstances(inst, view, count) {
  if (count <= 0) {
    return;
  }
  gl.bindBuffer(gl.ARRAY_BUFFER, inst.buffer);
  gl.bufferSubData(gl.ARRAY_BUFFER, 0, view, 0, count * inst.components);
  gl.bindVertexArray(inst.vao);
  gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, count);
  gl.bindVertexArray(null);
}

function drawGl(canvas) {
  gl.viewport(0, 0, canvas.width, canvas.height);
  let background = parseColor(c_background);
  gl.clearColor(background[0], background[1], background[2], 1);
  gl.clear(gl.COLOR_BUFFER_BIT);

  gl.useProgram(g_circle.prog);
  gl.uniform2f(g_circle.loc.u_resolution, canvas.width, canvas.height);

  { // Draw particles.
    g_particles = new Uint16Array(Module.HEAPU8.buffer, g_particles_ptr, g_particles_max_size);
    let size = GetParticles(g_particles.byteOffset, g_particles.length);
    gl.uniform1f(g_circle.loc.u_radius, 7);
    gl.uniform4fv(g_circle.loc.u_color, parseColor(c_green));
    drawInstances(g_gl_particles, g_particles, size >> 1);
  }

  gl.useProgram(g_line.prog);
  gl.uniform2f(g_line.loc.u_resolution, canvas.width, canvas.height);

  { // Draw portals, alternating colors between the two of each pair.
    g_portals = new Uint16Array(Module.HEAPU8.buffer, g_portals_ptr, g_portals_max_size);
    let size = GetPortals(g_portals.byteOffset, g_portals.length);
    gl.uniform1f(g_line.loc.u_width, 5);
    gl.uniform4fv(g_line.loc.u_color0, parseColor(c_blue));
    gl.uniform4fv(g_line.loc.u_color1, parseColor(c_orange));
    drawInstances(g_gl_portals, g_portals, size >> 2);
  }

  { // Draw bonds.
    g_bonds = new Uint16Array(Module.HEAPU8.buffer, g_bonds_ptr, g_bonds_max_size);
    let size = GetBonds(g_bonds.byteOffset, g_bonds.length);
    let color = parseColor("#ffffff80");
    gl.uniform1f(g_line.loc.u_width, 7);
    gl.uniform4fv(g_line.loc.u_color0, color);
    gl.uniform4fv(g_line.loc.u_color1, color);
    drawInstances(g_gl_bonds, g_bonds, size >> 2);
  }

  gl.useProgram(g_circle.prog);

  { // Draw frozen particles.
    g_frozen = new Uint16Array(Module.HEAPU8.buffer, g_frozen_ptr, g_frozen_max_size);
    let size = GetFrozen(g_frozen.byteOffset, g_frozen.length);
    gl.uniform1f(g_circle.loc.u_radius, 4);
    gl.uniform4fv(g_circle.loc.u_color, parseColor(c_black));
    drawInstances(g_gl_frozen, g_frozen, size >> 1);
  }
}

function draw() {
  let canvas = Module['canvas'];
  if (gl) {
    drawGl(canvas);
  } else {
    drawCanvas2d(canvas);
  }
}

// Fallback for browsers without WebGL2. Kept because a failed context would
// otherwise leave the page blank.
function drawCanvas2d(canvas) {
  let ctx = canvas.getContext('2d');
  ctx.drawImage(g_tmp_canvas, 0, 0, canvas.width, canvas.height);

  // Clear the canvas.
  ctx.fillStyle = c_background;
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  { // Draw particles.
    g_particles = new Uint16Array(Module.HEAPU8.buffer, g_particles_ptr, g_particles_max_size);
    let size = GetParticles(g_particles.byteOffset, g_particles.length);
    ctx.fillStyle = c_green;
    ctx.lineWidth = 0;
    radius = 7
    for (let i = 0; i + 1 < size; i += 2) {
      ctx.beginPath();
      ctx.arc(g_particles[i], g_particles[i + 1], radius, 0, 2 * Math.PI, true);
      ctx.fill();
    }
  }

  { // Draw portals.
    g_portals = new Uint16Array(Module.HEAPU8.buffer, g_portals_ptr, g_portals_max_size);
    let size = GetPortals(g_portals.byteOffset, g_portals.length);
    ctx.lineWidth = 5;
    for (let i = 0; i + 1 < size; i += 4) {
      xa = g_portals[i + 0];
      ya = g_portals[i + 1];
      xb = g_portals[i + 2];
      yb = g_portals[i + 3];
      // Do not draw empty portals, needed for incomplete pairs.
      if (xa != xb || ya != yb) {
        ctx.strokeStyle = (i / 4 % 2 == 0) ? c_blue : c_orange;
        ctx.beginPath();
        ctx.moveTo(xa, ya);
        ctx.lineTo(xb, yb);
        ctx.stroke();
      }
    }
  }

  { // Draw bonds.
    g_bonds = new Uint16Array(Module.HEAPU8.buffer, g_bonds_ptr, g_bonds_max_size);
    let size = GetBonds(g_bonds.byteOffset, g_bonds.length);
    ctx.lineWidth = 7;
    ctx.strokeStyle = "#ffffff80";
    for (let i = 0; i + 1 < size; i += 4) {
      ctx.beginPath();
      ctx.moveTo(g_bonds[i + 0], g_bonds[i + 1]);
      ctx.lineTo(g_bonds[i + 2], g_bonds[i + 3]);
      ctx.stroke();
    }
  }

  { // Draw frozen particles.
    g_frozen = new Uint16Array(Module.HEAPU8.buffer, g_frozen_ptr, g_frozen_max_size);
    let size = GetFrozen(g_frozen.byteOffset, g_frozen.length);
    ctx.fillStyle = c_black;
    ctx.lineWidth = 0;
    radius = 4
    for (let i = 0; i + 1 < size; i += 2) {
      ctx.beginPath();
      ctx.arc(g_frozen[i], g_frozen[i + 1], radius, 0, 2 * Math.PI, true);
      ctx.fill();
    }
  }
}

function restart() {
  Init();
  flag_pause = false;
  flag_accel = false;
  syncButtons();
}

function setButtonStyle(button_name, pressed) {
  document.getElementById('button_' + button_name).className = pressed ? "button pressed" : "button";
}

function syncButtons() {
  setButtonStyle('pause', flag_pause);
  setButtonStyle('accel', flag_accel);
  mousemode = GetMouseMode();
  setButtonStyle('r', mousemode == 'repulsion');
  setButtonStyle('a', mousemode == 'attraction');
  setButtonStyle('p', mousemode == 'pick');
  setButtonStyle('f', mousemode == 'freeze');
  setButtonStyle('o', mousemode == 'portal');
  setButtonStyle('b', mousemode == 'bonds');
  setButtonStyle('g', GetGravity());
}

function togglePause() {
  flag_pause = !flag_pause;
  syncButtons();
  SetPause(flag_pause);
}

function toggleAccel() {
  flag_accel = !flag_accel;
  if (flag_accel && !g_accel) {
    try {
      g_accel = new Accelerometer({ frequency: 5 });
      g_accel.addEventListener("reading", () => {
        if (flag_accel) {
          SetGravityVect(-g_accel.x, -g_accel.y);
          gx = -g_accel.x.toFixed(3);
          gy = -g_accel.y.toFixed(3);
          window.text_accel.innerHTML = `<br>gravity=(${gx}, ${gy})`;
        }
      });
      g_accel.start();
    } catch (error) {
      console.log(error);
    }
  }
  syncButtons();
}


function clearOutput() {
  if (output) {
    output.value = '';
  }
  if (outputerr) {
    outputerr.value = '';
  }
}
function print(text) {
  if (output) {
    output.value += text + "\n";
    output.scrollTop = output.scrollHeight;
  }
}
function printError(text) {
  if (outputerr) {
    outputerr.value += text + "\n";
    outputerr.scrollTop = outputerr.scrollHeight;
  }
}

// Codes for keys without one, matching Control::kKeyArrow* in control.h.
const kKeyArrowDown = 256;
const kKeyArrowUp = 257;

function sendKeyDownChar(c) {
  if (c == 'ArrowDown') {
    SendKeyDown(kKeyArrowDown);
    return;
  }
  if (c == 'ArrowUp') {
    SendKeyDown(kKeyArrowUp);
    return;
  }
  keysym = c.charCodeAt(0);
  if (c.length == 1 && keysym < 256) {
    SendKeyDown(keysym);
  }
}

function pressButton(c) {
  if (c == 'q') {
    restart();
  } else {
    sendKeyDownChar(c);
    syncButtons();
  }
}

function postRun() {
  GetParticles = Module.cwrap('GetParticles', 'number', ['number', 'number']);
  GetPortals = Module.cwrap('GetPortals', 'number', ['number', 'number']);
  GetBonds = Module.cwrap('GetBonds', 'number', ['number', 'number']);
  GetFrozen = Module.cwrap('GetFrozen', 'number', ['number', 'number']);

  SendKeyDown = Module.cwrap('SendKeyDown', null, ['number']);
  SendMouseMotion = Module.cwrap('SendMouseMotion', null, ['number', 'number']);
  SendMouseDown = Module.cwrap('SendMouseDown', null, ['number', 'number']);
  SendMouseUp = Module.cwrap('SendMouseUp', null, ['number', 'number']);
  SetControlDebug = Module.cwrap('SetControlDebug', null, ['number']);
  SetPause = Module.cwrap('SetPause', null, ['number']);
  GetGravity = Module.cwrap('GetGravity', 'number', []);
  SetGravity = Module.cwrap('SetGravity', null, ['number']);
  SetGravityVect = Module.cwrap('SetGravityVect', null, ['number', 'number']);
  GetMouseMode = Module.cwrap('GetMouseMode', 'string', []);
  Init = Module.cwrap('Init', null, []);

  g_particles_ptr = Module._malloc(g_particles_max_size * 2);
  g_portals_ptr = Module._malloc(g_portals_max_size * 2);
  g_bonds_ptr = Module._malloc(g_bonds_max_size * 2);
  g_frozen_ptr = Module._malloc(g_frozen_max_size * 2);

  let canvas = Module['canvas'];
  // A canvas keeps the first context type it is asked for, so this decides
  // once which of the two renderers draw() uses.
  if (!initGl(canvas)) {
    console.log('WebGL2 unavailable, falling back to the 2D canvas');
    g_tmp_canvas = document.createElement('canvas');
    g_tmp_canvas.width = canvas.width;
    g_tmp_canvas.height = canvas.height;
  }

  let handler_keydown = function(e) {
    if (e.key == ' ') {
      togglePause();
      return;
    }
    pressButton(e.key);
    // Prevent scrolling by Space and by the arrow keys.
    if((e.keyCode == 32 || e.key == 'ArrowDown' || e.key == 'ArrowUp') &&
       e.target == document.body) {
      e.preventDefault();
    }
  };
  let get_xy = function(e) {
    let x = -1 + 2 * e.offsetX / canvas.clientWidth;
    let y = 1 - 2 * e.offsetY / canvas.clientHeight;
    return [x, y]
  };
  let get_touch_xy = function(e, touch) {
    var rect = e.target.getBoundingClientRect();
    var x = -1 + 2 * (e.changedTouches[0].pageX - rect.left) / canvas.clientWidth;
    var y = 1 - 2 * (e.changedTouches[0].pageY - rect.top) / canvas.clientHeight;
    return [x, y]
  };


  let handler_mousemove = function(e) {
    e.preventDefault();
    xy = get_xy(e);
    SendMouseMotion(xy[0], xy[1]);
  };
  let handler_mousedown = function(e) {
    e.preventDefault();
    xy = get_xy(e);
    SendMouseDown(xy[0], xy[1]);
  };
  let handler_mouseup = function(e) {
    e.preventDefault();
    xy = get_xy(e);
    SendMouseUp(xy[0], xy[1]);
  };
  let handler_touchmove = function(e) {
    e.preventDefault();
    xy = get_touch_xy(e);
    SendMouseMotion(xy[0], xy[1]);
  };
  let handler_touchstart = function(e) {
    e.preventDefault();
    xy = get_touch_xy(e);
    SendMouseDown(xy[0], xy[1]);
  };
  let handler_touchend = function(e) {
    e.preventDefault();
    xy = get_touch_xy(e);
    SendMouseUp(xy[0], xy[1]);
  };

  window.addEventListener('keydown', handler_keydown);
  canvas.addEventListener('mousemove', handler_mousemove);
  canvas.addEventListener('mousedown', handler_mousedown);
  canvas.addEventListener('mouseup', handler_mouseup);
  canvas.addEventListener('touchmove', handler_touchmove);
  canvas.addEventListener('touchstart', handler_touchstart);
  canvas.addEventListener('touchend', handler_touchend);

  // Disable Space on buttons.
  [
    window.button_pause, window.button_restart,
    window.button_r, window.button_a, window.button_p,
    window.button_f, window.button_o, window.button_b,
    window.button_i, window.button_g, window.button_d,
  ].forEach(b => {
    b.addEventListener('keydown', function(e){
      if (e.key == ' ') {
        e.preventDefault();
      }
    }, false);
    b.addEventListener('keyup', function(e){
      if (e.key == ' ') {
        e.preventDefault();
      }
    }, false);
  });

  syncButtons();
}

var Module = {
  preRun: [],
  postRun: [postRun],
  print: (function(text) {
    clearOutput();
    return function(text) {
      if (arguments.length > 1) {
        text = Array.prototype.slice.call(arguments).join(' ');
      }
      print(text);
    };
  })(),
  printErr: (function(text) {
    clearOutput();
    return function(text) {
      if (arguments.length > 1) {
        text = Array.prototype.slice.call(arguments).join(' ');
      }
      printError(text);
    };
  })(),
  canvas: (function() { return document.getElementById('canvas'); })(),
  setStatus: function(text) {},
};
