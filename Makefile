all: \
	index.html \

PANDOC = pandoc
PFLAGS =
STYLE = .

.md.html:
	$(PANDOC) -s --css $(STYLE)/pandoc.css --katex $(PFLAGS) "$<" -o "$@"

clean:
	rm -vf *.pdf *.html

.PHONY: all clean
.SUFFIXES: .md .html
