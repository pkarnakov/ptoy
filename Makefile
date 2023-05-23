BUILDDIR = build
MAKEFILE = $(BUILDDIR)/Makefile
CMAKE = cmake

all: $(MAKEFILE)
	+make -C $(BUILDDIR) $(target)

$(MAKEFILE):
	mkdir -p "$(BUILDDIR)"
	(cd "$(BUILDDIR)" && $(CMAKE) -DUSE_AVX=1 ..)

web:
	git clone -b gh-pages --single-branch git@github.com:pkarnakov/ptoy.git web

clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean install
