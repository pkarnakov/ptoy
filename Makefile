BUILD = build
MAKEFILE = $(BUILD)/Makefile
CMAKE = cmake

WASM_BUILD = build_wasm
WASM_MAKEFILE = $(WASM_BUILD)/Makefile
WASM_CMAKE = emcmake cmake
WASM_FLAGS = \
 -DUSE_OPENMP=0 -DUSE_BACKEND_SDL=0 -DUSE_AVX=0 \
 -DUSE_WASM=1 -DUSE_EXE=0 -DUSE_EXEM=0

all: $(MAKEFILE)
	+make -C $(BUILD) $(target)

$(MAKEFILE):
	mkdir -p "$(BUILD)"
	(cd "$(BUILD)" && $(CMAKE) -DUSE_AVX=1 ..)

$(WASM_MAKEFILE): Makefile
	mkdir -p "$(WASM_BUILD)"
	cd "$(WASM_BUILD)" && $(WASM_CMAKE) .. $(WASM_FLAGS)

wasm: $(WASM_MAKEFILE)
	+make -C $(WASM_BUILD) $(target)

serve:
	cd "$(WASM_BUILD)" && emrun --serve_after_exit ptoy.html

pages:
	git clone -b gh-pages --single-branch git@github.com:pkarnakov/ptoy.git pages

clean:
	rm -rf $(BUILD) $(WASM_BUILD)

.PHONY: all wasm serve clean
