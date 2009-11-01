#################################################################
# Infernal
# src Makefile
# SVN $Id: Makefile.in 2560 2008-08-25 14:30:12Z nawrockie $
#################################################################
# @LICENSE@
#################################################################

SHELL  = /bin/sh

## your compiler and compiler flags
#
CC     = gcc -std=gnu99
CFLAGS = -g -O0 -mpentiumpro -Wall
SIMDFLAGS= -msse2

## other defined flags. 
#  DEFS contains stuff that autoconf decides on. 
#  LIBS contains system libraries that the configure script decides we need.
#
DEFS  = -DHAVE_CONFIG_H
LIBS  =  -lm

# Where the build arena is (module .h and .o files)
#
MYLIBS   = -linfernal -lhmmer -leasel
MYLIBDIR = -L../ -L../../easel -L../../hmmer/src/ 
MYINCDIR = -I../ -I../../easel -I../../hmmer/src/

PROGS =

OBJS  = cm_optimized.o\
        sse_cm_dpsearch.o\
        sse_cm_dpsmall.o\
        sse_cmcons_hitmx.o\
        sse_cmcons_mscyk.o\
        sse_util.o

HDRS  = ../funcs.h\
	../structs.h\
        impl_sse.h

.c.o: 
	$(CC) $(SIMDFLAGS) $(CFLAGS) $(DEFS) $(MYINCDIR) -c $<

#################################################################
## Targets defining how to make INFERNAL executables.
##
all: 	$(OBJS) 

$(PROGS): %: %.o $(OBJS) ${HDRS}
	$(CC) $(SIMDFLAGS) $(CFLAGS) $(DEFS) $(MYINCDIR) $(MYLIBDIR) -o $@ $@.o $(OBJS) $(MYLIBS) $(LIBS)

test: $(OBJS) ${HDRS}
	$(CC) $(SIMDFLAGS) $(CFLAGS) $(DEFS) $(MYINCDIR) $(MYLIBDIR) -o sse-cyk -DIMPLSSE_SMALL_TEST sse_cm_dpsmall.c sse_util.o cm_optimized.o $(MYLIBS) $(LIBS)
	$(CC) $(SIMDFLAGS) $(CFLAGS) $(DEFS) $(MYINCDIR) $(MYLIBDIR) -o sse-bmark -DIMPL_SEARCH_BENCHMARK sse_cm_dpsearch.c $(MYLIBS) $(LIBS)
	$(CC) $(SIMDFLAGS) $(CFLAGS) $(DEFS) $(MYINCDIR) $(MYLIBDIR) -o sse-mscyk -DIMPL_MSCYK_TEST sse_cmcons_mscyk.c sse_cmcons_hitmx.o sse_cm_dpsearch.o cm_optimized.o $(MYLIBS) $(LIBS)
	$(CC) $(SIMDFLAGS) $(CFLAGS) $(DEFS) $(MYINCDIR) $(MYLIBDIR) -o mscyk-test -DIMPL_MSCYK_TEST2 sse_cmcons_mscyk.c sse_cmcons_hitmx.o sse_cm_dpsearch.o cm_optimized.o $(MYLIBS) $(LIBS)
	

#################################################################
## Miscellaneous targets.
##
clean:
	-rm -f *.o *~ Makefile.bak core $(PROGS) TAGS gmon.out libinfernal.a
	-rm -f sse-bmark sse-cyk sse-mscyk

distclean:
	-rm -f *.o *~ Makefile.bak core $(PROGS) TAGS gmon.out libinfernal.a
	-rm -f Makefile config.h

TAGS:
	etags -t *.c *.h Makefile.in

