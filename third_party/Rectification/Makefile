# Generated automatically from Makefile.in by configure.
CXX = c++
CCMALLOC = ccmalloc
CFLAGS =  -O3 -DNDEBUG -pipe
CPPFLAGS = -O3 -DNDEBUG
CXXFLAGS = -Wall
LDFLAGS = -ggdb -lm -L/afs/cs/u/aliensun/scratch/view3/third_party/CLAPACK -L/afs/cs/u/aliensun/scratch/view3/third_party/CLAPACK/F2CLIBS -llapack_LINUX -lblas_LINUX -lF77
%LDFLAGS = -ggdb -lm -llapack -lblas -lf2c
srcdir = .
OBJECTS = src/main.o src/Rectify/Ops.o src/Rectify/General2Im.o

RUNFLAGS:= Test/Image0.pnm Test/Image1.pnm Test/FM Test/matches Test/RectIm1.pnm Test/RectIm2.pnm

test: Rectification
	./Rectification ${RUNFLAGS}

Rectification: ${OBJECTS}
	${CXX} ${OBJECTS} ${LDFLAGS} -o Rectification

efence: ${OBJECTS}
	${CXX} ${OBJECTS} ${LDFLAGS} -lefence -o Rectification
	echo "set args Rectification ${RUNFLAGS}" > gdbcommands
	ddd Rectification -x gdbcommands

.PHONY : clean
clean:
	@rm -f $(srcdir)/src/*.o $(srcdir)/src/*~ $(srcdir)/src/*% ${OBJECTS} TAGS
	@rm -f $(srcdir)/src/Rectify/*.o $(srcdir)/src/Rectify/*~ $(srcdir)/src/Rectify/*% ${OBJECTS} TAGS

.PHONY : distclean
distclean:
	@make clean
	@cd $(srcdir); rm -f *.a *.bak *~ *% Rectification gdbcommands; rm -Rf auto*cache TAGS BROWSE Test/RectIm1.pnm Test/RectIm2.pnm VisStudioDotNetProject/debug  VisStudioDotNetProject/release
	$(srcdir)/cleanup .

.PHONE : package
package:
	@make distclean
	@rm -f ../Rectification.tar.bz2
	@cd ..; tar -cvf Rectification.tar Rectification; bzip2 -9 Rectification.tar

.PHONY : debug
debug: ${OBJECTS}
	echo "set args ${RUNFLAGS}" > gdbcommands
	ddd Rectification -x gdbcommands

include Makefile.depend
