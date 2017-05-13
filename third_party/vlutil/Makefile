# file:        Makefile
# author:      Andrea Vedaldi
# description: Build SIFT mex files

DIST=vlutil
VER=1.1.0
BINDIST=$(DIST)-$(VER)-$(ARCH)

# --------------------------------------------------------------------
#
# --------------------------------------------------------------------

# Determine on the flight the system we are running on
Darwin_ARCH := mac
Linux_ARCH  := glx
ARCH := $($(shell uname)_ARCH)

#mac_CFLAGS       := -O -I.  -pedantic -Wall
mac_CFLAGS       := -O -I.  
mac_MEX_CFLAGS   := CFLAGS='$$CFLAGS $(mac_CFLAGS)'
mac_MEX_SUFFIX   := mexmac

# glx_CFLAGS       := -O -I. -pedantic -Wall
glx_CFLAGS       := -O -I. 
glx_MEX_CFLAGS   := CFLAGS='$$CFLAGS $(glx_CFLAGS)'
glx_MEX_SUFFIX   := mexglx

MEX_SUFFIX   := $($(ARCH)_MEX_SUFFIX)
MEX_CFLAGS   := $($(ARCH)_MEX_CFLAGS)

# --------------------------------------------------------------------
#
# --------------------------------------------------------------------

vpath %.mex.c   toolbox
vpath %.mex.cpp toolbox

src  := $(wildcard toolbox/*.mex.c) $(wildcard toolbox/*.mex.cpp)
msrc := $(wildcard toolbox/*.m)
stem := $(notdir $(basename $(basename $(src))))
tgt  := $(addprefix mex/, $(addsuffix .$(MEX_SUFFIX),$(stem)))

mex/%.$(MEX_SUFFIX) : %.mex.c
	mex -I. $(MEX_CFLAGS) $< -outdir 'mex'
	@mv mex/$*.mex.$(MEX_SUFFIX) mex/$*.$(MEX_SUFFIX)

mex/%.$(MEX_SUFFIX) : %.mex.cpp
	mex -I. $(MEX_CFLAGS) $< -outdir 'mex'
	@mv mex/$*.mex.$(MEX_SUFFIX) mex/$*.$(MEX_SUFFIX)

.PHONY: all
all: $(tgt)

.PHONY: info
info :
	@echo src  = $(src)
	@echo stem = $(stem)
	@echo tgt  = $(tgt)

# PDF documentation
.PHONY: doc
doc: vlutil.html

vlutil.html : $(msrc)
	mdoc --output=vlutil.html --exclude 'test_*' toolbox

.PHONY: clean
clean:
	rm -f $(tgt)
	rm -f .DS_Store .gdb_history
	find . -name '*~' -exec rm -f \{\} \;

.PHONY: distclean
distclean: clean
	rm -f mex/*.mexmac mex/*.mexglx
	rm -f vlutil.html
	rm -rf $(DIST)-*

.PHONY: dist
dist: distclean
	echo Version $(VER) >TIMESTAMP
	echo Archive created on `date` >>TIMESTAMP
	tar chzvf $(DIST)-$(VER).tar.gz ../$(notdir $(CURDIR))

.PHONY: bindist
bindist: all doc
	test -e $(BINDIST) || mkdir $(BINDIST)
	cp README toolbox/*.m vlutil.html mex/*.$(MEX_SUFFIX) $(BINDIST)
	cd $(BINDIST) ; strip -S *.$(MEX_SUFFIX) $(notdir $(LIB))	
	echo Version $(VER) >$(BINDIST)/TIMESTAMP
	echo Package created on `date` >>$(BINDIST)/TIMESTAMP
	tar cvzf $(BINDIST).tar.gz $(BINDIST)

.PHONY: autorights
autorights: distclean
	autorights \
	  . \
	  --recursive \
	  --template gpl \
	  --years 2006 \
	  --authors "Andrea Vedaldi" \
	  --holders "Andrea Vedaldi" \
	  --program "VLUtil"
