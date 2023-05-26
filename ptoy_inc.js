var SetConfig;
var GetConfig;
var output = document.getElementById('output');
var outputerr = document.getElementById('outputerr');
var g_tmp_canvas;
var g_particles_max_size = 10000;
var kScale = 2;

var GetParticles;
var g_particles;
var g_particles_ptr;

function Draw() {
  if (output) {
    output.value = GetConfig();
  }

  let canvas = Module['canvas'];
  let ctx = canvas.getContext('2d');
  ctx.drawImage(g_tmp_canvas, 0, 0, canvas.width, canvas.height);

  ctx.clearRect(0, 0, canvas.width, canvas.height);

  ctx.lineWidth = 3;
  ctx.strokeStyle="#000000";
  ctx.strokeRect(0, 0, canvas.width, canvas.height);

  g_particles = new Uint16Array(Module.HEAPU8.buffer, g_particles_ptr, g_particles_max_size);
  let size = GetParticles(g_particles.byteOffset, g_particles.length);
  ctx.fillStyle = "#00CD6C";
  ctx.lineWidth = 0;
  radius = 3.5
  for (let i = 0; i + 1 < size; i += 2) {
    ctx.beginPath();
    ctx.arc(g_particles[i], g_particles[i + 1], radius, 0, 2 * Math.PI, true);
    ctx.fill();
  }
}

function ClearOutput() {
  output.value = '';
  outputerr.value = '';
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
