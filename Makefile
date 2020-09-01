BUILDDIR = build
MAKEFILE = $(BUILDDIR)/Makefile
CMAKE = cmake

all: $(MAKEFILE)
	+make -C $(BUILDDIR) $(target)

$(MAKEFILE):
	mkdir -p "$(BUILDDIR)"
	(cd "$(BUILDDIR)" && $(CMAKE) -DUSE_AVX=1 ..)

clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean install
