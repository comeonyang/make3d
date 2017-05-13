
# ------------------ Compilation options ------------------------

# Loads math library.  
LIBS = -lm
GET = get
#CFLAGS = -Wall -march=Xeon -mfpmath=sse -mmmx -msse -msse2 -msse3 -ansi
CFLAGS = -Wall -march=pentium4 -mfpmath=sse -mmmx -msse -msse2 -msse3 -ansi
#CFLAGS = -Wall -O3 -DNO_DEBUG -march=pentium4 -mfpmath=sse -mmmx -msse -msse2 -msse3 -ansi
DEPEND= makedepend $(CFLAGS)

CC = g++
CXX = g++

#CC = g++ -fabi-version=4.0.2
#CXX = g++ -fabi-version=4.0.2

# --------------------- Code modules ----------------------------

# Source files
SRCS = main.cpp imload.cpp os_mapping.cpp
# Object files
OBJ = main.o imload.o os_mapping.o
# Definitions
DEFS = image.h fasthessian.h ipoint.h surf.h imload.h

# ------------------------ Rules --------------------------------
#$(SRCS):
#        $(GET) $@

# Link against static library
surf.ln: ${OBJ} libSurf.so
	${CC} -o $@ ${CFLAGS} main.o imload.o os_mapping.o -static libSurf.a ${LIBS}

# Small matching demo application

match.ln: match.cpp imload.o 
	${CC} -o $@ ${CFLAGS} imload.o match.cpp -static libSurf.a -lm
symmatch.ln: symmatch.cpp imload.o 
	${CC} -o $@ ${CFLAGS} imload.o symmatch.cpp -static libSurf.a -lm
symmatchConS.ln: symmatchConS.cpp imload.o 
	${CC} -o $@ ${CFLAGS} imload.o symmatchConS.cpp -static libSurf.a -lm
symmatchRConS.ln: symmatchRConS.cpp imload.o 
	${CC} -o $@ ${CFLAGS} imload.o symmatchRConS.cpp -static libSurf.a -lm
symmatchRConS_thre.ln: symmatchRConS_thre.cpp imload.o 
	${CC} -o $@ ${CFLAGS} imload.o symmatchRConS_thre.cpp -static libSurf.a -lm
symDenseMatch.ln: symDenseMatch.cpp imload.o 
	${CC} -o $@ ${CFLAGS} imload.o symDenseMatch.cpp -static libSurf.a -lm
OcluMatchRConS_thre.ln: OcluMatchRConS_thre.cpp imload.o
	${CC} -o $@ ${CFLAGS} imload.o OcluMatchRConS_thre.cpp -static libSurf.a -lm

# To link against a shared library, use
#surf.ln: ${OBJ} libSurf.so
#	${CC} -o $@ ${CFLAGS} main.o imload.o -L. -lSurf ${LIBS}
# Note to set LD_LIBRARY_PATH environment variable before running surf.ln
all: surf.ln match.ln symmatch.ln symmatchConS.ln symmatchRConS.ln symDenseMatch.ln symmatchRConS_thre.ln OcluMatchRConS_thre.ln

clean:
	-rm *.o surf.ln match.ln symmatch.ln symDenseMatch.ln symmatchConS.ln symmatchRConS.ln symmatchRConS_thre.ln OcluMatchRConS_thre.ln

#depend: $(SRCS)
#        $(DEPEND) $(SRCS)
