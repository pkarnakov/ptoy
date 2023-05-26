var SetConfig;
var GetConfig;
var output = document.getElementById('output');
var outputerr = document.getElementById('outputerr');
var g_tmp_canvas;
var kScale = 2;

var GetParticles;
var g_particles;
var g_particles_ptr;
var g_particles_max_size = 10000;

var GetPortals;
var g_portals;
var g_portals_ptr;
var g_portals_max_size = 10000;

var c_greed = "#00cd6c";
var c_blue = "#009ade";
var c_orange = "#f28522";

function Draw() {
  let canvas = Module['canvas'];
  let ctx = canvas.getContext('2d');
  ctx.drawImage(g_tmp_canvas, 0, 0, canvas.width, canvas.height);

  // Clear the canvas.
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  // Draw particles.
  {
    g_particles = new Uint16Array(Module.HEAPU8.buffer, g_particles_ptr, g_particles_max_size);
    let size = GetParticles(g_particles.byteOffset, g_particles.length);
    ctx.fillStyle = c_greed;
    ctx.lineWidth = 0;
    radius = 7
    for (let i = 0; i + 1 < size; i += 2) {
      ctx.beginPath();
      ctx.arc(g_particles[i], g_particles[i + 1], radius, 0, 2 * Math.PI, true);
      ctx.fill();
    }
  }

  // Draw portals.
  {
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
}

function ClearOutput() {
  if (output) {
    output.value = '';
  }
  if (outputerr) {
    outputerr.value = '';
  }
}
function PrintError(text) {
  console.error(text);
  if (outputerr) {
    outputerr.value += text + "\n";
    outputerr.scrollTop = outputerr.scrollHeight;
  }
}
function PostRun() {
  SetConfig = Module.cwrap('SetConfig', 'int', ['string']);
  GetConfig = Module.cwrap('GetConfig', 'string', []);
  GetParticles = Module.cwrap('GetParticles', 'number', ['number', 'number']);
  GetPortals = Module.cwrap('GetPortals', 'number', ['number', 'number']);

  g_particles_ptr = Module._malloc(g_particles_max_size * 2);

  let canvas = Module['canvas'];
  g_tmp_canvas = document.createElement('canvas');
  g_tmp_canvas.width = canvas.width;
  g_tmp_canvas.height = canvas.height;
}

var Module = {
  preRun: [],
  postRun: [PostRun],
  printErr: (function(text) {
    ClearOutput();
    return function(text) {
      if (arguments.length > 1) {
        text = Array.prototype.slice.call(arguments).join(' ');
      }
      PrintError(text);
    };
  })(),
  canvas: (function() { return document.getElementById('canvas'); })(),
  setStatus: function(text) {},
};
