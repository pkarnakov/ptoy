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


// Misc
var SetPause;
var g_pause = false;
var GetGravity;
var SetGravity;
var GetMouseMode;
var Init;


function draw() {
  let canvas = Module['canvas'];
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
    for (let i = 0; i + 1 < size; i += 8) {
      ctx.strokeStyle = c_blue;
      ctx.beginPath();
      ctx.moveTo(g_portals[i + 0], g_portals[i + 1]);
      ctx.lineTo(g_portals[i + 2], g_portals[i + 3]);
      ctx.stroke();

      ctx.strokeStyle = c_orange;
      ctx.beginPath();
      ctx.moveTo(g_portals[i + 4], g_portals[i + 5]);
      ctx.lineTo(g_portals[i + 6], g_portals[i + 7]);
      ctx.stroke();
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
    g_frozen = new Uint16Array(Module.HEAPU8.buffer, g_bonds_ptr, g_bonds_max_size);
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
  g_pause = false;
  syncButtons();
}

function setButtonStyle(button_name, pressed) {
  document.getElementById('button_' + button_name).className = pressed ? "button pressed" : "button";
}

function syncButtons() {
  setButtonStyle('pause', g_pause);
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
  g_pause = !g_pause;
  syncButtons();
  SetPause(g_pause);
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

function sendKeyDownChar(c) {
  keysym = c.charCodeAt(0);
  if (keysym < 256) {
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
  GetMouseMode = Module.cwrap('GetMouseMode', 'string', []);
  Init = Module.cwrap('Init', null, []);

  g_particles_ptr = Module._malloc(g_particles_max_size * 2);
  g_portals_ptr = Module._malloc(g_portals_max_size * 2);
  g_bonds_ptr = Module._malloc(g_bonds_max_size * 2);
  g_frozen_ptr = Module._malloc(g_frozen_max_size * 2);

  let canvas = Module['canvas'];
  g_tmp_canvas = document.createElement('canvas');
  g_tmp_canvas.width = canvas.width;
  g_tmp_canvas.height = canvas.height;

  let handler_keydown = function(e) {
    if (e.key == ' ') {
      togglePause();
      return;
    }
    pressButton(e.key);
    // Prevent scrolling by Space.
    if(e.keyCode == 32 && e.target == document.body) {
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
    window.button_i, window.button_g,
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
