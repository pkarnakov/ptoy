# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

ptoy is a 2D interactive particle simulation (Lennard-Jones-like pair forces, bonds,
portals, mouse forces). It builds as a native SDL2/OpenGL app and as a WebAssembly
page rendered by JavaScript on a 2D canvas.

## Build and run

```
make                 # Release build into build/ (USE_AVX=1), binary build/ptoy
make debug           # Debug build into build_debug/, same options, assertions on
make target=<t>      # forward a target to the generated build/Makefile (also with debug)
make wasm            # emcmake build into build_wasm/ (no OpenMP, no SDL; AVX via SIMD128)
make serve           # emrun build_wasm/ptoy.html
make clean           # remove build/, build_debug/ and build_wasm/
./build/ptoy         # runnable from any directory, see PTOY_ASSETS_DIR below
```

Use `make debug` when a crash or a memory error is suspected: it keeps `assert`
active, notably the alignment and capacity checks around the AVX kernel.

There is no test suite and no linter config beyond `.clang-format` (Google-based,
80 columns, 2-space indent). `USE_SANDBOX=1` builds `src/sandbox.cpp`, a scratch
executable for experiments.

CMake options that select code paths: `USE_AVX`, `USE_OPENMP`, `USE_EXE`,
`USE_WASM`, `USE_BACKEND_SDL`, `USE_BACKEND_TEXT`, `USE_SANDBOX`,
`USE_WARNINGS`, `USE_TRAPEZOID`. Reconfigure an existing build directory with
`(cd build && cmake -DUSE_TRAPEZOID=1 ..)` and then `make` as usual.
In C++ they are tested with `USEFLAG(X)` from `src/macros.h`, which expands to
`0USE_X` so an undefined flag evaluates to 0 — write `#if USEFLAG(BACKEND_SDL)`,
not `#ifdef USE_BACKEND_SDL`.

## Architecture

Two entry points share the same simulation core:

- `src/main.cpp` — native. Owns the frame loop, calls `View::Control()` then
  `display()`, which advances the simulation and copies state into a `Scene`.
- `src/ptoy_wasm.cpp` — WebAssembly. Same loop under `emscripten_set_main_loop`,
  but instead of drawing it exposes `extern "C"` getters (`GetParticles`,
  `GetPortals`, `GetBonds`, `GetFrozen`) that pack pixel coordinates into
  caller-provided `uint16_t` buffers, plus input setters (`SendKeyDown`,
  `SendMouse*`, `SetGravity*`, `SetPause`). `src/ptoy_inc.js` `cwrap`s these,
  reads the buffers as `Uint16Array` views on `Module.HEAPU8`, and draws them
  with WebGL2. Any new export must be added to both the `EXPORTED_FUNCTIONS`
  list in `CMakeLists.txt` and the `cwrap` block in `ptoy_inc.js`.

Layers:

- `Particles` (`src/particles.{h,cpp}`) — the simulation. Velocity-Verlet-style
  half steps in `Particles::step()`, parallelized with OpenMP over blocks.
  Also owns portals, bonds, frozen particles, and the mouse force; the mouse
  interactions are exposed as `XxxStart/XxxMove/XxxStop(Vect)` triples.
- `blocks` (`src/blocks.h`) — uniform grid of cells (`kBlockSize = 4 * kRadius`)
  holding struct-of-arrays particle data in `AlignedAllocator` vectors.
  `SortParticles()` re-bins particles after each step and drops those that left
  the domain. `block_by_id_` maps particle id to (block, index).
- `Control` (`src/control.{h,cpp}`) — keyboard/mouse events to `Particles` calls,
  holds `MouseMode`. Shared by the SDL and wasm frontends; the key bindings in
  `Control::Handle` are the single source of truth for both (see README for the list).
- `View` / `Scene` (`src/scene.h`) — the rendering interface. `Scene` is a bundle
  of non-owning `span`s; the frontend fills `SceneData` vectors each frame and
  hands out spans into them. Implementations: `ViewGl` (`src/view_gl.cpp`,
  pimpl over SDL2 + OpenGL, also draws the GUI buttons and text) and `ViewText`
  (`src/view_text.h`, prints to stdout, used with `USE_BACKEND_TEXT`).
- `Game` (`src/game.h`) — owns the `Particles` instance and maps window size to
  the simulation domain (`[-1,-1]` to `[1,1]` at the reference 800x800 size).

`Scal` is `float` and `Vect` is a 2-component struct (`src/geometry.h`); the whole
simulation is single precision.

## Things that bite

- **AVX kernel.** `CalcForceAvx` in `src/particles.cpp` processes 8 particles at a
  time with aligned loads/stores and runs past `size()` up to the next multiple of
  8. It works because block arrays use `AlignedAllocator<Vect, 64>` and are
  pre-`reserve`d (`kBlockPadding` in `blocks.h`). Changing the block data
  containers, the allocator, or the reserve size can silently break it. The
  `CALC_FORCE` macro selects `CalcForceAvx` or `CalcForceSerial` at compile time.
- **Shader and asset paths** go through `Asset()` in `view_gl.cpp`, which prefixes
  them with `PTOY_ASSETS_DIR` — the source directory recorded at build time by
  CMake — or with `$PTOY_ASSETS` if set. The executable therefore runs from any
  working directory, but it points at the tree it was built from.
- **The WebGL renderer in `ptoy_inc.js`** binds the packed `uint16_t` buffers
  straight as instanced vertex attributes (`gl.UNSIGNED_SHORT`, unnormalized),
  so the layout written by `GetParticles` and friends is part of its interface —
  changing the packing or the entry size breaks drawing with no error. Circles
  and lines are both instanced unit quads; lines are expanded in the vertex
  shader because WebGL caps `gl_LineWidth` at 1. A canvas keeps the first
  context type requested, so `initGl()` decides once between WebGL2 and the
  `drawCanvas2d` fallback; the two must stay visually in sync.
- **Font assets** (`assets/font/font.{bin,geom,png}`) are generated by
  `assets/font/font.py` (PIL + numpy) and read back by `LoadFont` in `view_gl.cpp`;
  regenerate both files together.
- **Physical constants** (`kRadius`, `kSigma*`, `kMass`, `dt`, `kVelocityLimit`)
  live at the top of `src/particles.cpp` and are tightly coupled — `kBlockSize`
  must stay at least the force cutoff, since forces are only computed between
  neighboring blocks.
- **The integrator in `Particles::step()`** is a midpoint force evaluation with
  a first-order position update by default; `-DUSE_TRAPEZOID=1` selects the
  trapezoidal one, second order in both variables, and `t` toggles it at
  runtime. `calc_forces()` takes no velocity parameter: the corrector needs the
  dashpot evaluated at $v^\ast$, and gets it only because the predictor
  overwrote `data.velocity` in place first. Reordering those loops, or using a
  saved copy of the velocity, silently drops the scheme to first order. See
  `docs/MODEL.md`.
- **Contact damping comes from the dashpot in `F12()`** (`kDashpot`), not from
  `kDissipation` and no longer from integrator error. It acts only along the
  line of centres and only while particles overlap, so it settles piles without
  slowing bulk flow. Removing it makes piles shiver indefinitely.
- **Assertions** use `fassert`, `fassert_equal`, `NAMEVALUE` from `src/logger.h`
  (they throw, they are not `assert`).

## Publishing the web version

`pages/` is a separate clone of this repo on the `gh-pages` branch
(`make pages` creates it). From `pages/`: `make wasm` rsyncs the artifacts out of
`build_wasm/` and records the source revision, then `make commit`.
