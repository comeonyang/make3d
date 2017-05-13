/* src/config.h.  Generated automatically by configure.  */
// Size of certain variables

// relational operators are screwed - so make sure don't get defined and mess everything up in older versions of gcc (pre 3.0).
#define __SGI_STL_INTERNAL_RELOPS

#define SIZEOF_UNSIGNED_CHAR 1
#define SIZEOF_CHAR 1
#define SIZEOF_UNSIGNED_SHORT 2
#define SIZEOF_SHORT 2
#define SIZEOF_UNSIGNED_LONG_INT 4
#define SIZEOF_LONG_INT 4

// Can link with the bsd library for GNU C.
// Allows us to switch on GNU extensions and get loads of juicy
// inline optimisations and extra optimised functions
#define HAVE_LIBBSD 1
#ifdef __GNUC__
#if HAVE_LIBBSD==1
// On better set up systems _GNU_SOURCE will already be defined for me
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif
#endif
