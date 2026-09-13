BUILD = build
MAKEFILE = $(BUILD)/Makefile
CMAKE = cmake

DEBUG_BUILD = build_debug
DEBUG_MAKEFILE = $(DEBUG_BUILD)/Makefile

WASM_BUILD = build_wasm
WASM_MAKEFILE = $(WASM_BUILD)/Makefile
WASM_CMAKE = emcmake cmake
WASM_FLAGS = \
 -DUSE_OPENMP=0 -DUSE_BACKEND_SDL=0 -DUSE_AVX=0 -DUSE_WASM=1 -DUSE_EXE=0

default: cmake

cmake: $(MAKEFILE)
	+make -C $(BUILD) $(target)

$(MAKEFILE):
	mkdir -p "$(BUILD)"
	(cd "$(BUILD)" && $(CMAKE) -DUSE_AVX=1 -DCMAKE_BUILD_TYPE=Release ..)

debug: $(DEBUG_MAKEFILE)
	+make -C $(DEBUG_BUILD) $(target)

$(DEBUG_MAKEFILE):
	mkdir -p "$(DEBUG_BUILD)"
	(cd "$(DEBUG_BUILD)" && $(CMAKE) -DUSE_AVX=1 -DCMAKE_BUILD_TYPE=Debug ..)

$(WASM_MAKEFILE):
	mkdir -p "$(WASM_BUILD)"
	cd "$(WASM_BUILD)" && $(WASM_CMAKE) .. $(WASM_FLAGS)

wasm: $(WASM_MAKEFILE)
	+make -C $(WASM_BUILD) $(target)

serve:
	cd "$(WASM_BUILD)" && emrun --serve_after_exit ptoy.html

pages:
	git clone -b gh-pages --single-branch git@github.com:pkarnakov/ptoy.git pages

clean:
	rm -rf $(BUILD) $(DEBUG_BUILD) $(WASM_BUILD)

.PHONY: default cmake debug wasm serve clean
