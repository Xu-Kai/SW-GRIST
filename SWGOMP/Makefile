all: lib/libswgomp.a plugin/swgomp.so
.PHONY: lib/libswgomp.a  plugin/swgomp.so clean clean-lib clean-plugin doc/_build/latex/swgomp-doc.pdf
PREFIX ?= ./
lib/libswgomp.a:
	+make -C lib
plugin/swgomp.so:
	+make -C plugin
docpdf: doc/_build/latex/swgomp-doc.pdf
	cp $< swgomp-doc.pdf
doc/_build/latex/swgomp-doc.pdf:
	+make -C doc latexpdfja
clean: clean-lib clean-plugin
clean-lib:
	+make -C lib clean
clean-plugin:
	+make -C plugin clean
PREFIX ?= ./
SCRIPTS := $(wildcard script/*)
MODULES := $(wildcard lib/fortran/*.mod)
H_CROSS := $(wildcard include/cross/*.h)
H_SLAVE := $(wildcard include/slave/*.h)
install: lib/libswgomp.a plugin/swgomp.so $(SCRIPTS) $(MODULES) $(H_CROSS) $(H_SLAVE)
	mkdir -p $(PREFIX)/bin/
	mkdir -p $(PREFIX)/lib/
	mkdir -p $(PREFIX)/plugin/
	mkdir -p $(PREFIX)/include/cross/
	mkdir -p $(PREFIX)/include/slave/
	mkdir -p $(PREFIX)/include/fortran/
	cp -P $(SCRIPTS) $(PREFIX)/bin/
	cp plugin/swgomp.so $(PREFIX)/plugin/
	cp lib/libswgomp.a $(PREFIX)/lib/
	cp $(H_CROSS) $(PREFIX)/include/cross/
	cp $(H_SLAVE) $(PREFIX)/include/slave/
	cp $(MODULES) $(PREFIX)/include/fortran/
