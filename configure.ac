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
AC_INIT(Infernal, 1.0rc1, eddys@janelia.hhmi.org, infernal)
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

PACKAGE_DATE="June 2008"
PACKAGE_COPYRIGHT="Copyright (C) 2008 HHMI Janelia Farm Research Campus"
PACKAGE_LICENSE="Freely distributed under the GNU General Public License (GPL)"
PACKAGE_LICENSETAG=gnu
PACKAGE_VERSION=$PACKAGE_VERSION
PACKAGE_URL="http://infernal.janelia.org/"

EASEL_DATE="June 2008"
EASEL_COPYRIGHT="Copyright (C) 2008 HHMI Janelia Farm Research Campus"
EASEL_LICENSE="Freely distributed under the Janelia Software License."
EASEL_LICENSETAG=jsl
EASEL_VERSION="i1.0rc1"
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

AC_ARG_ENABLE(gcov,[AS_HELP_STRING([--enable-gcov],[compile for code coverage testing])], enable_gcov=$enableval,  enable_gcov=no)
AC_ARG_ENABLE(sse,[AS_HELP_STRING([--enable-sse],[enable SSE optimizations])] ,           enable_sse=$enableval,   enable_sse=check)
AC_ARG_ENABLE(vmx,[AS_HELP_STRING([--enable-vmx],[enable Altivec/VMX optimizations])],    enable_vmx=$enableval,   enable_vmx=check)
AC_ARG_WITH(gsl,[AS_HELP_STRING([--with-gsl],[use the GSL, GNU Scientific Library])],     with_gsl=$withval,       with_gsl=no)
AC_ARG_ENABLE(mpi,[AS_HELP_STRING([--enable-mpi],[enable MPI parallelization])],          enable_mpi=$enableval,   enable_mpi=no)

# end of optionally configured features
################################################################

################################################################
# 4. Checks for programs
################################################################

AC_PROG_CC
AC_PROG_CPP
AX_COMPILER_VENDOR
#  AC_PROG_CC_STDC wants to append -c99 to icc; icc doesn't like this.
if test "$ax_cv_c_compiler_vendor" != "intel"; then
  AC_PROG_CC_STDC
fi
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])
CHECK_GNU_MAKE

# Select our default optimization flags in CFLAGS.
#
if test "$sre_cflags_env_set" = "no"; then
if test "$enable_gcov" = "yes"; then
   # --enable-gcov overrides CC too.
   CC="gcc"
   CFLAGS="-g -Wall -fprofile-arcs -ftest-coverage"
elif test "$enable_debugging" != "no"; then
   if test "$GCC" = "yes"; then
      CFLAGS="-g -Wall"
   fi
else
   sre_save_ac_test_cflags="$ac_test_CFLAGS"
   if test "$sre_cflags_env_set" = "no"; then
     ac_test_CFLAGS=""
   fi
   AX_CC_MAXOPT
   ac_test_CFLAGS="$sre_save_ac_test_cflags"
   # That was a workaround. AX_CC_MAXOPT apparently thinks that
   # $ac_test_CFLAGS only gets set if ${CFLAGS} was set by the
   # user (in the environment), but that's not so; it gets set
   # by AC_PROG_CC. Instead of rewriting AX_CC_MAXOPT, we fool
   # it.
fi
fi

# MPI :  set @CC@ to mpicc, sets @MPILIBS@, defines HAVE_MPI. 
if test "$enable_mpi" = "yes"; then
   ACX_MPI(,AC_MSG_ERROR([MPI library not found for --enable-mpi]))
   CC=$MPICC
fi

# SSE: figure out our @SIMD_CFLAGS@, define HAVE_SSE2
if test "$enable_sse" = "check" && test "$enable_vmx" != "yes"; then
   AC_CHECK_HEADER([emmintrin.h],[enable_sse=yes])
fi
if test "$enable_sse" = "yes"; then
   if test "x$SIMD_CFLAGS" = x; then
     case $ax_cv_c_compiler_vendor in
       gnu)   AX_CHECK_COMPILER_FLAGS(-msse2,[SIMD_CFLAGS="-msse2"],[]);;
       *)     ;;
     esac   
   fi
   AC_DEFINE(HAVE_SSE2,1,[Support SSE2 (Streaming SIMD Extensions 2) instructions])
fi

# VMX/Altivec (not autodetected yet; must use --enable-vmx to enable)
if test "$enable_vmx" = "yes"; then
   if test "x$SIMD_CFLAGS" = x; then
   if test "$GCC" = "yes"; then
     AX_CHECK_COMPILER_FLAGS(-faltivec, [SIMD_CFLAGS="-faltivec"],
       [AX_CHECK_COMPILER_FLAGS(-maltivec -mabi=altivec,
 	   [SIMD_CFLAGS="-maltivec -mabi=altivec"],
           [AX_CHECK_COMPILER_FLAGS(-fvec, [SIMD_CFLAGS="-fvec"],
		   [AC_MSG_ERROR([Need a version of gcc with -maltivec])])])])
   elif test "${ax_cv_c_compiler_vendor}" = "ibm"; then
     AX_CHECK_COMPILER_FLAGS(-qaltivec, [SIMD_CFLAGS="-qaltivec"],[])
   else
     AX_CHECK_COMPILER_FLAGS(-faltivec, [SIMD_CFLAGS="-faltivec"],[])
   fi
   AC_DEFINE(HAVE_VMX,1,[Support Altivec/VMX vector instructions])
   fi
fi
AC_SUBST(SIMD_CFLAGS)


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
AC_CHECK_HEADERS([sys/types.h unistd.h stdint.h inttypes.h])

# altivec.h requires the simd cflags
# For reasons I don't understand, this needs to come after any other CHECK_HEADERS().
if test "$enable_vmx" = "yes"; then
   sre_save_CFLAGS="$CFLAGS"
   sre_save_CPPFLAGS="$CPPFLAGS"
   CFLAGS="$CFLAGS $SIMD_CFLAGS"
   CPPFLAGS="$CPPFLAGS $SIMD_CFLAGS"
   AC_CHECK_HEADERS([altivec.h])
   CFLAGS="$sre_save_CFLAGS"
   CPPFLAGS="$sre_save_CPPFLAGS"
fi

# Checks for types
#
AC_C_BIGENDIAN
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
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

