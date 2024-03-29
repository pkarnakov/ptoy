# ptoy

[Play web version](https://pkarnakov.github.io/ptoy/ptoy.html)

| Native | [Web](https://pkarnakov.github.io/ptoy/ptoy.html) |
|:---:|:---:|
|<img src="https://pkarnakov.github.io/ptoy/images/ptoy_native.png" height="170">|<img src="https://pkarnakov.github.io/ptoy/images/ptoy_web.png" height="170"> |


Screenshots from older version
|   |   |   |
|:---:|:---:|:---:|
|<img src="https://pkarnakov.github.io/ptoy/images/screenshot0.png" height="100">|<img src="https://pkarnakov.github.io/ptoy/images/screenshot1.png" height="100">|<img src="https://pkarnakov.github.io/ptoy/images/screenshot2.png" height="100">|

## Clone

```
git clone https://github.com/pkarnakov/ptoy.git
```

## Requirements

* Linux or Mac OS X
* C++14 compiler
* CMake (optional)
* emscripten (optional to build WebAssembly)
* SDL2
* glew

## Build

Build with CMake

```
make
```

Build without CMake, using `pkgconf` to find SDL2 and glew

```
make legacy
```

## Run

```
cd build && ./ptoy
```

## Control

### Keyboard

* `r`: switch mouse to *Repulsion* mode
* `a`: switch mouse to *Attraction* mode
* `p`: switch mouse to *Pick* mode
* `f`: switch mouse to *Freeze* mode
* `o`: switch mouse to *Portal* mode
* `b`: switch mouse to *Bonds* mode
* `n`: switch mouse to *No action* mode
* `g`: toggle gravity
* `i`: remove last pair of portals
* `space`: toggle pause
* `q`: quit after pressing 3 times

### Mouse

Click actions depend on the mode selected by keyboard,
all mouse buttons are equivalent.

Mode|Key|Click action
:---|:---:|---
*Repulsion* (default) |`r`| repulsive force
*Attraction* |`a`| attractive force
*Portal*  |`o`| draw portals (blue, then orange)
*Pick*  |`p`| pick the closest particle
*Freeze* |`f`| freeze/unfreeze particle
*Bonds*  |`b`| draw bonds between particles
*No action* |`n`| nothing
