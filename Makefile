all: \
	info.html \

PANDOC = pandoc
PFLAGS =
STYLE = .
REPO = ..
RSYNC = rsync -a -i --update

.md.html:
	$(PANDOC) -s --css $(STYLE)/pandoc.css --katex $(PFLAGS) "$<" -o "$@"

wasm:
	(cd "$(REPO)" && git rev-parse --short HEAD) > revision
	$(RSYNC) $(REPO)/build_wasm/{ptoy{.html,.js,.wasm,_inc.js,.css},favicon.png,libs} .

commit:
	git commit -m "update from devel `cat revision`"

clean:
	rm -vf *.pdf *.html

.PHONY: all clean wasm
.SUFFIXES: .md .html
