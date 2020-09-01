BUILDDIR = build
MAKEFILE = $(BUILDDIR)/Makefile
CMAKE = cmake

all: $(MAKEFILE)
	+make -C $(BUILDDIR) $(target)

$(MAKEFILE):
	mkdir -p "$(BUILDDIR)"
	(cd "$(BUILDDIR)" && $(CMAKE) ..)

clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean install
