# Process this file with autoconf to produce the HMMER configure script.
#
# Note that this is synchronized with Easel's configure script,
# so that HMMER configures Easel in addition to itself; this saves
# having to do a separate ./configure in Easel. That is,
# everything that appears in Easel's configure script must also
# appear here. 
#
# reminders to save re-reading autoconf manual for the n'th time:
#   - output variables:
#        are defined here as normal shell variables, e.g. FOO="my string"
#        are made into output variables by calling AC_SUBST(FOO)
#        any occurrence of @FOO@ in an output file is then substituted
#        I think this only happens in files we assign w/ AC_CONFIG_FILES;
#        that's the Makefile.in's.
#
#   - C preprocessor symbols:
#        are defined here by calling AC_DEFINE(FOO) or AC_DEFINE(FOO, [42])
#        then #undef FOO lines in easel.h.in become #define FOO or #define FOO 42
#        I think this only happens in header files that we assign
#        w/ AC_CONFIG_HEADERS -- which means, easel.h.in
#
# SRE, Mon Oct  5 14:55:45 1998
# SVN $Id$
# xref autoconf macro archive: //www.gnu.org/software/ac-archive/
#
# GNU recommends the following order:
#   1. autoconf requirements
#   2. AC_INIT
#   3. info on the package
#   4. checks for programs
#   5. checks for libraries
#   6. checks for header files
#   7. checks for types
#   8. checks for structures
#   9. checks for compiler characteristics
#  10. checks for library functions
#  11. checks for system services
#  12. AC_CONFIG_FILES
#  13. AC_OUTPUT


AC_PREREQ(2.60)
AC_INIT(Infernal, 1.0, eddys@janelia.hhmi.org, infernal)
AC_MSG_NOTICE([Configuring Infernal for your system.])


################################################################
# 3. Info on the package
#
# The four AC_INIT args set these output variables and preprocessor symbols:
#     PACKAGE_NAME      <package>     e.g. "Infernal"
#     PACKAGE_VERSION   <version>     e.g. "1.0"
#     PACKAGE_BUGREPORT <bug-report>  e.g. "eddys@janelia.hhmi.org"
#     PACKAGE_TARNAME   <tarname>     e.g. "infernal"
# From them, it derives one more:
#     PACKAGE_STRING    <package> <version>, e.g. "Infernal 1.0"
# We also define additional variables:
#     PACKAGE_RELCODE     unique CVS tag without .'s: e.g. "infernal1_0"
#     PACKAGE_DATE        release date: e.g. "October 2007"
#     PACKAGE_COPYRIGHT   one-line copyright string
#     PACKAGE_LICENSE     one-line license string
#     PACKAGE_LICENSETAG  which license to bundle from Licenses/ subdirectory.
#     PACKAGE_VERSION     copy of version code, e.g. "1.0"
#     PACKAGE_URL         URL home.
# And we have to define the relevant package variables for Easel as well.
#
################################################################

PACKAGE_RELCODE="infernal1_0"
PACKAGE_DATE="October 2007"
PACKAGE_COPYRIGHT="Copyright (C) 2001-2007 HHMI Janelia Farm Research Campus"
PACKAGE_LICENSE="Freely distributed under the GNU General Public License (GPL)"
PACKAGE_LICENSETAG=gnu
PACKAGE_VERSION=$PACKAGE_VERSION
PACKAGE_URL="http://infernal.janelia.org/"

EASEL_RELCODE="easel0_1"
EASEL_DATE="June 2007"
EASEL_COPYRIGHT="Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus"
EASEL_LICENSE="Freely distributed under the Janelia Software License."
EASEL_LICENSETAG=jsl
EASEL_VERSION="0.1"
EASEL_URL="http://easel.janelia.org/"

# Output variables  (AC_OUTPUT replaces @var@ in input files, such as Makefiles)
AC_SUBST(PACKAGE_RELCODE)
AC_SUBST(PACKAGE_DATE)
AC_SUBST(PACKAGE_COPYRIGHT)
AC_SUBST(PACKAGE_LICENSE)
AC_SUBST(PACKAGE_LICENSETAG)
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(PACKAGE_URL)

AC_SUBST(EASEL_RELCODE)
AC_SUBST(EASEL_DATE)
AC_SUBST(EASEL_COPYRIGHT)
AC_SUBST(EASEL_LICENSE)
AC_SUBST(EASEL_LICENSETAG)
AC_SUBST(EASEL_VERSION)
AC_SUBST(EASEL_URL)

# Preprocessor symbols (replace #undefs in p7config.h)
AC_DEFINE_UNQUOTED(PACKAGE_DATE,      "$PACKAGE_DATE")
AC_DEFINE_UNQUOTED(PACKAGE_COPYRIGHT, "$PACKAGE_COPYRIGHT")
AC_DEFINE_UNQUOTED(PACKAGE_LICENSE,   "$PACKAGE_LICENSE")
AC_DEFINE_UNQUOTED(PACKAGE_VERSION,   "$PACKAGE_VERSION")

AC_DEFINE_UNQUOTED(EASEL_DATE,      "$EASEL_DATE")
AC_DEFINE_UNQUOTED(EASEL_COPYRIGHT, "$EASEL_COPYRIGHT")
AC_DEFINE_UNQUOTED(EASEL_LICENSE,   "$EASEL_LICENSE")
AC_DEFINE_UNQUOTED(EASEL_VERSION,   "$EASEL_VERSION")

AC_DEFINE(eslLIBRARY)

# Figure out what host we're compiling on.
# Three GNU scripts must be included in the distro: 
#       install.sh, config.guess, config.sub
# This sets four shell variables:
#       host            example: i686-pc-linux-gnu      
#       host_cpu        example: i686
#       host_vendor     example: pc
#       host_os         example: linux-gnu
#
# These are used later in the configure by the ACX_PTHREAD macro to determine
# system-specific details of threads libraries.
#
AC_CANONICAL_HOST

# Check if we're a binary distribution. 
# Trigger is existence of binaries/, nonexistence of src/.
# 
# For binary distro, we do minimal configuration, setting version information 
# and installation path names, making ./configure independent of having a working
# C compiler.
#
# Implements a sneaky if statement that wraps most of the configure.ac M4 code.
#
if ((! test -d binaries) && test -d src); then 
AC_MSG_NOTICE([Full Infernal source distribution - using full configuration])

# Choose the compiler, and get it in ISO ANSI C (C99) mode.
#
AC_PROG_CC
if test "$CC" = xlc ; then	# IBM's xlc?
  CFLAGS=$CFLAGS" -O3 -qarch=auto -qtune=auto"
fi
# Work around a known bug in autoconf's AC_PROG_CC_C99's handling of icc
if test "$CC" != icc; then 
  AC_PROG_CC_C99
else
  CFLAGS=$CFLAGS" -std=c99"
fi

################################################################
# Process the ./configure command line
################################################################

# --enable-debugging=x    - set debugging level to <x> (1-3)
#
# At all levels, including 0, replaces CFLAGS w/ "-g -Wall" (so it assumes gcc).
# Sets the p7_DEBUGLEVEL preprocessor symbol to <x>
#
AC_ARG_ENABLE(debugging,
[  --enable-debugging      set CFLAGS for debugging
  --enable-debugging=x    also set diagnostics level to <x> (1-3) (3 = most verbose)],
[ case $enable_debugging in
   yes)  AC_MSG_NOTICE([enabled debugging diagnostics level 0 (CFLAGS only, no verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(cm_DEBUGLEVEL, 0)
	 AC_DEFINE(eslDEBUGLEVEL, 0)
         ;;
     1)  AC_MSG_NOTICE([enabled debugging diagnostics level 1 (low verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(cm_DEBUGLEVEL, 1)
	 AC_DEFINE(eslDEBUGLEVEL, 1)
         ;;
     2)  AC_MSG_NOTICE([enabled debugging diagnostics level 2 (moderate verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(cm_DEBUGLEVEL, 2)
	 AC_DEFINE(eslDEBUGLEVEL, 2)
         ;;
     3)  AC_MSG_NOTICE([enabled debugging diagnostics level 3 (high verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(cm_DEBUGLEVEL, 3)
	 AC_DEFINE(eslDEBUGLEVEL, 3)
         ;;
    no)  AC_MSG_NOTICE([debugging diagnostics disabled])
         AC_DEFINE(cm_DEBUGLEVEL, 0)
	 AC_DEFINE(eslDEBUGLEVEL, 0)
         ;;
     *)  echo "Ignoring unknown argument to --enable-debugging: $enable_debugging"
         ;;
esac])

# --enable-gcov      - compile for code coverage testing
#
# Replaces CC with "gcc" and CFLAGS with "-g -Wall -fprofile-arcs -ftest-coverage".
# Running programs (notably the test suite) will then generate .gcda files,
# which the gcov tool can read.
#
AC_ARG_ENABLE(gcov, 
  [  --enable-gcov           Compile with instrumentation for code coverage testing], 
  [ 
    AC_MSG_NOTICE([Compiling with gcov instrumentation.])
    CC="gcc"
    CFLAGS="-g -Wall -fprofile-arcs -ftest-coverage"
  ])

# --enable-mpi       Enable MPI parallelization
# 
# Sets MPILIBS output variable
# Sets CC to appropriate compiler wrapper for MPI (mpicc, usually)
# Defines HAVE_MPI preprocessor variable.
#
AC_ARG_ENABLE(mpi, 
[AS_HELP_STRING([--enable-mpi], [enable MPI parallelization])],
[case $enable_mpi in
   yes) ACX_MPI([
           AC_MSG_NOTICE([enabled optional MPI parallelization])
	   AC_DEFINE(HAVE_MPI)
	   CC=$MPICC
          ], [
           AC_MSG_ERROR([MPI library not found])
          ])
	 ;;
   no)   AC_MSG_NOTICE([MPI parallelization disabled])
         ;;
   *)    echo "Ignoring unknown argument to --enable-mpi: $enable_mpi"
	 ;;
esac])

# --enable-sse      Force SSE2 SIMD optimization on
# --disable-sse     Force SSE2 SIMD optimization off
#                   Default is to have SSE2 on, if it's detected.
#
# This sets @SIMD_FLAGS@ output variable, defines HAVE_SSE2, and
# if SSE is chosen, sets impl_choice to "sse".
AC_ARG_ENABLE(sse,
   [AS_HELP_STRING([--enable-sse],  [ enable SSE SIMD parallelization (default)]) ],
   [
     if test "x$withval" != xno ; then
        # icc is going to generate SSE code anyway, and it'll complain about -msse2 flag
       if test "$CC" != icc; then 
         AC_DEFINE(HAVE_SSE2, 1)
         AX_EXT
       else 
         impl_choice="sse"
       fi
     else
       AC_DEFINE(HAVE_SSE2, 0)
     fi
   ],
   [
     AX_EXT
   ])
if test "$ax_have_sse2_ext" = yes; then
  impl_choice="sse"
else
  AC_MSG_NOTICE([no SSE detected])
fi

# --with-gsl        - enable hooks into the GSL (GNU Scientific Library)
#
# The (action-if-given) is blank: using --with-gsl sets $with_gsl to 'yes',
# and --without-gsl sets $with_gsl to 'no'. If neither is given,
# $with_gsl is set to 'check', and we'll try to use GSL anyway if we
# can find it.
# 
# We interpret the $with_gsl setting later, in the "checks for libraries" section.
AC_ARG_WITH([gsl],
            [AS_HELP_STRING([--with-gsl],
                           [use the GSL, GNU Scientific Library (default is no)])],
            [if test "x$withval" != xno ; then 
               AC_MSG_NOTICE([Enabling hooks into the GNU Scientific Library (GSL)])
             fi
            ],
	    [with_gsl=no])

# --enable-lfs           Large File Summit (LFS) support for >2GB files
# See: http://ftp.sas.com/standards/large.file/x_open.20Mar96.html
#
AC_ARG_ENABLE(lfs, 
[  --enable-lfs            enable LFS, Large File Support],
[case $enable_lfs in
   yes) AC_MSG_NOTICE([configured for optional LFS, large file support])
	AC_DEFINE(_LARGEFILE_SOURCE)
	AC_DEFINE(_LARGEFILE64_SOURCE)
	AC_DEFINE(_FILE_OFFSET_BITS, 64)
        ;;
   no)  ;;
    *)  echo "Ignoring unknown argument to --enable-lfs: $enable_lfs"
        ;;
esac])

## Rigoroous filters not included as an option in v1.0+
## Users can use previous versions of Infernal inside which
## rigfilters did successfully build.
##
## Default option for rigfilters is OFF
#rigfilters=0;
## --with-rigfilters     install rigorous filters
#AC_ARG_WITH(rigfilters,
#[  --with-rigfilters       enable rigorous filters, requires CFSQP package],
#[case $with_rigfilters in
#   yes) AC_MSG_NOTICE([Rigorous filters selected])
#	rigfilters=1
#	;;
#   no)  rigfilters=0
#        ;;
#    *)  echo "Ignoring unknown argument to --with-rigfilters: $with_rigfilters"
#        ;;
#esac])
#
## Default location of cfsqp
#cfsqp_path=$PWD/rigfilters/cfsqp;
## --with-cfsqp             specify location of cfsqp
#AC_ARG_WITH(cfsqp,
#[  --with-cfsqp=PATH       specify the location of CFSQP source files],
#[case $with_cfsqp in
#   /*)  cfsqp_path=$with_cfsqp
#        ;;
#   ~*)  cfsqp_path=`echo $with_cfsqp | sed 's/^~//'`
#        cfsqp_path=$HOME/../$cfsqp_path
#        ;;
#   *)   #NOTE: relative path not normalized
#        cfsqp_path=$PWD/$with_cfsqp
#        ;;
#  esac
##  if test -d $cfsqp_path; then
##    AC_MSG_NOTICE([path recognized as directory])
##  else
##    AC_MSG_NOTICE([path not recognized as directory!: $cfsqp_path])
##  fi
#  if test -f $cfsqp_path/cfsqp.c; then
#    cp $cfsqp_path/cfsqp.c rigfilters/cfsqp/.
#  fi
#  if test -f $cfsqp_path/cfsqpusr.h; then
#    cp $cfsqp_path/cfsqpusr.h rigfilters/cfsqp/.
#  fi
#  if test -f $cfsqp_path/qld.c; then
#    cp $cfsqp_path/qld.c rigfilters/cfsqp/.
#  fi
#])
#if test "$rigfilters" = 1; then
#  if test -f rigfilters/cfsqp/cfsqp.c; then
#    :
#  else
#    AC_MSG_NOTICE([CFSQP: Unable to locate cfsqp.c])
#    rigfilters=0;
#  fi
#  if test -f rigfilters/cfsqp/cfsqpusr.h; then
#    :
#  else
#    AC_MSG_NOTICE([CFSQP: Unable to locate cfsqpusr.h])
#    rigfilters=0
#  fi
#  if test -f rigfilters/cfsqp/qld.c; then
#    :
#  else
#    AC_MSG_NOTICE([CFSQP: Unable to locate qld.c])
#    rigfilters=0
#  fi
#  if test "$rigfilters" = 0; then
#    AC_MSG_NOTICE([CFSQP: Source not found; disabling rigorous filters])
#    rigfilters=""
#  else
#    AC_MSG_NOTICE([CFSQP located, including rigfilters package in build])
#    rigfilters="rigfilters"
#  fi
#fi
#AC_SUBST(RIGFILTERS,$rigfilters)
#
## Default option for rigfilters is OFF
#rigfilters=0;
## --with-rigfilters     install rigorous filters
#AC_ARG_WITH(rigfilters,
#[  --with-rigfilters       enable rigorous filters, requires CFSQP package],
#[case $with_rigfilters in
#   yes) AC_MSG_NOTICE([Rigorous filters selected])
#	rigfilters=1
#	;;
#   no)  rigfilters=0
#        ;;
#    *)  echo "Ignoring unknown argument to --with-rigfilters: $with_rigfilters"
#        ;;
#esac])
#
# end of rigfilters block
#################################

# end of optionally configured features
################################################################

# Set output variable "RANLIB"
AC_PROG_RANLIB

# Set output variable "AR"
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])

# Set the output variable EXEC_DEPENDENCY according
# to whether we're using GNU make or a SYSV make.
CHECK_GNU_MAKE



# Checks for libraries
#
LIBGSL=
AS_IF([test "x$with_gsl" != xno],
      [AC_CHECK_LIB([gsl], [gsl_expm1], 
           [AC_SUBST([LIBGSL], ["-lgsl -lgslcblas"])
            AC_DEFINE([HAVE_LIBGSL], [1], [Define if you have libgsl])
           ],
           [if test "x$with_gsl" != xcheck; then
             AC_MSG_FAILURE(
               [--with-gsl was given, but GSL library was not found])
            fi
           ],
           [-lgslcblas]
        )])

# Checks for headers
#
AC_CHECK_HEADERS([sys/types.h])
AC_CHECK_HEADERS([unistd.h])
AC_CHECK_HEADERS([stdint.h])
AC_CHECK_HEADERS([inttypes.h])

# Checks for types
#
AC_C_BIGENDIAN
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_OFF_T

# Checks for functions, defining HAVE_FOO when foo is found
#
AC_CHECK_FUNCS(mkstemp)
AC_CHECK_FUNCS(popen)
AC_CHECK_FUNCS(strcasecmp)
AC_CHECK_FUNCS(times)
AC_FUNC_FSEEKO
AC_CHECK_FUNCS(ntohs, , AC_CHECK_LIB(socket, ntohs))
AC_CHECK_FUNCS(ntohl, , AC_CHECK_LIB(socket, ntohl))
AC_CHECK_FUNCS(htons, , AC_CHECK_LIB(socket, htons))
AC_CHECK_FUNCS(htonl, , AC_CHECK_LIB(socket, htonl))

#
# 11. System services
#
AC_SYS_LARGEFILE

# Done.
# Config subdirs and files (except main Makefile, which we defer)
#

AC_CONFIG_HEADERS(src/config.h)

# Write Infernal's config.h header
#
AC_CONFIG_FILES(src/Makefile)
AC_CONFIG_FILES(testsuite/Makefile)
AC_CONFIG_FILES(rigfilters/Makefile)
AC_CONFIG_FILES(rigfilters/cfsqp/Makefile)
AC_CONFIG_FILES(rigfilters/cm2hmm-1.0/Makefile)

AC_CONFIG_HEADERS([easel/esl_config.h])
AC_CONFIG_FILES([easel/miniapps/Makefile])
AC_CONFIG_FILES([easel/testsuite/Makefile])
AC_CONFIG_FILES([easel/Makefile])


################################################################
# 13. AC_OUTPUT
################################################################
# remember that sneaky if statement that wraps the source
# configuration? well, it ends here.
else
  AC_MSG_NOTICE([This is a precompiled binary distribution - using abbreviated config])
fi
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

