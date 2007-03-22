# Process this file with autoconf to produce a configure script.
# SRE, Thu Jul 27 13:30:41 2000
# SVN $Id$
# Hint: See autoconf macro archive at 
#     http://www.gnu.org/software/ac-archive/

# Initialize w/ AC_INIT(<package> <version> <bug-report> <tarname>)
# 
# The four AC_INIT args set the following output variables and preprocessor symbols:
#     PACKAGE_NAME      <package>     e.g. "Infernal"
#     PACKAGE_VERSION   <version>     e.g. "0.42"
#     PACKAGE_BUGREPORT <bug-report>  e.g. "eddy@genetics.wustl.edu"
#     PACKAGE_TARNAME   <tarname>     e.g. "infernal"
# From them, it derives one more:
#     PACKAGE_STRING    <package> <version>, e.g. "Infernal 0.42"
# We also define additional variables:
#     PACKAGE_RELCODE     unique CVS tag without .'s: e.g. "infernal0_42"
#     PACKAGE_DATE        release date: e.g. "September 2003"
#     PACKAGE_COPYRIGHT   one-line copyright string
#     PACKAGE_LICENSE     one-line license string
#     PACKAGE_LICENSETAG  which license to bundle from Licenses/ subdirectory.
#
# The preprocessor symbols get picked up in config.h.
# The output variables are used in Makefiles.
#
AC_PREREQ(2.57)
AC_INIT(Infernal, 0.72, eddys@janelia.hhmi.org, infernal)
AC_MSG_NOTICE([Welcome. Configuring Infernal for your system.])

PACKAGE_RELCODE="infernal0_72"
PACKAGE_DATE="January 2007"
PACKAGE_COPYRIGHT="Copyright (C) 2001-2007 HHMI Janelia Farm"
PACKAGE_LICENSE="Freely distributed under the GNU General Public License (GPL)"
PACKAGE_LICENSETAG=gnu

AC_SUBST(PACKAGE_RELCODE)
AC_SUBST(PACKAGE_DATE)
AC_SUBST(PACKAGE_COPYRIGHT)
AC_SUBST(PACKAGE_LICENSE)
AC_SUBST(PACKAGE_LICENSETAG)

AC_DEFINE_UNQUOTED(PACKAGE_DATE, "$PACKAGE_DATE")
AC_DEFINE_UNQUOTED(PACKAGE_COPYRIGHT, "$PACKAGE_COPYRIGHT")
AC_DEFINE_UNQUOTED(PACKAGE_LICENSE, "$PACKAGE_LICENSE")

# Check if we're a full source distribution, as opposed to a binary distro.
# The trigger for this is the nonexistence of the binaries/ subdirectory,
# and the existence of a src/ subdirectory.
# 
# If we are a binary dist, we only do a minimal configuration, setting up
# version information and installation path names.
# This makes ./configure independent of having a working
# C compiler.
#
# Watch out: this is implemented as a sneaky if statement that wraps the
# majority of the configure.ac M4 code.
#
if ((! test -d binaries) && test -d src); then 
AC_MSG_NOTICE([This is a full source distribution - using full config])


# Checks for programs.
# Sets output variables CC, CFLAGS, and RANLIB.
# AR might be hidden away on Solaris; provide some hints in the PATH.
#
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])

# GNU make check; makefile dependency lines for executables.
# original from John Darrington <j.darrington@elvis.murdoch.edu.au> 
# w/ heavy modifications from me.
# We need this because GNU make and SYSV make use different systems
# specifying variables for dependencies: $$@ in sysv, %: %.o in GNU.
# Would love to hear a better way of doing this.
# 
AC_DEFUN(CHECK_GNU_MAKE,[ 
  AC_MSG_CHECKING(whether your make is GNU make)
  foundGNUmake='nope.' ;
  EXEC_DEPENDENCY=[\$\$\@.o] ;
  if ( make --version nothing 2> /dev/null | grep GNU > /dev/null ) ;  then
     foundGNUmake='yes, it is.' ;
     EXEC_DEPENDENCY='%: %.o' ;
  fi
  AC_MSG_RESULT($foundGNUmake)
  AC_SUBST(EXEC_DEPENDENCY)
])
CHECK_GNU_MAKE

################################################################
# Now we're into our optional configuration flags
#


# --enable-ccmalloc      - configure for compiling in ccmalloc memory checks
#
# The CC output variable is modified.
# Requires ccmalloc version >= 0.4.0 because of the --no-wrapper option.
#
AC_ARG_ENABLE(ccmalloc,
[  --enable-ccmalloc       turn ccmalloc memory debugging on ],
[case $enable_ccmalloc in
   yes)  AC_MSG_NOTICE([enabled ccmalloc memory debugging])
         CC="ccmalloc --no-wrapper $CC" 
         ;;
    no)  AC_MSG_NOTICE([ccmalloc memory debug/tracing disabled])
         ;;
     *)  echo "Ignoring unknown argument to --enable-ccmalloc: $enable_ccmalloc"
	 ;;
esac])



# --enable-debugging=x    - set debugging level to <x> (1-3)
# 
# At all levels, including 0, replaces CFLAGS w/ "-g -Wall" (so it assumes gcc).
# Sets the DEBUGLEVEL preprocessor symbol to <x>
#
AC_ARG_ENABLE(debugging,
[  --enable-debugging      set CFLAGS for debugging
  --enable-debugging=x    also set diagnostics level to <x> (1-3) 
                          (where 3 = most verbose)],
[ case $enable_debugging in
   yes)  AC_MSG_NOTICE([enabled debugging diagnostics level 0 (CFLAGS only, no verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(DEBUGLEVEL, 0)
         ;;
     1)  AC_MSG_NOTICE([enabled debugging diagnostics level 1 (low verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(DEBUGLEVEL, 1)
         ;;
     2)  AC_MSG_NOTICE([enabled debugging diagnostics level 2 (moderate verbosity)])
         CFLAGS="-g -Wall"
	 AC_DEFINE(DEBUGLEVEL, 2)
         ;;
     3)  AC_MSG_NOTICE([enabled debugging diagnostics level 3 (high verbosity)])
         CFLAGS="-g -Wall"
	 AC_DEFINE(DEBUGLEVEL, 3)
         ;;
    no)  AC_MSG_NOTICE([debugging diagnostics disabled])
         ;;
     *)  echo "Ignoring unknown argument to --enable-debugging: $enable_debugging"
	 ;;
esac])

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

# Default option for rigfilters is OFF
rigfilters=0;
# --with-rigfilters     install rigorous filters
AC_ARG_WITH(rigfilters,
[  --with-rigfilters       enable rigorous filters, requires CFSQP package],
[case $with_rigfilters in
   yes) AC_MSG_NOTICE([Rigorous filters selected])
	rigfilters=1
	;;
   no)  rigfilters=0
        ;;
    *)  echo "Ignoring unknown argument to --with-rigfilters: $with_rigfilters"
        ;;
esac])

# Default location of cfsqp
cfsqp_path=$PWD/rigfilters/cfsqp;
# --with-cfsqp             specify location of cfsqp
AC_ARG_WITH(cfsqp,
[  --with-cfsqp=PATH       specify the location of CFSQP source files],
[case $with_cfsqp in
   /*)  cfsqp_path=$with_cfsqp
        ;;
   ~*)  cfsqp_path=`echo $with_cfsqp | sed 's/^~//'`
        cfsqp_path=$HOME/../$cfsqp_path
        ;;
   *)   #NOTE: relative path not normalized
        cfsqp_path=$PWD/$with_cfsqp
        ;;
  esac
#  if test -d $cfsqp_path; then
#    AC_MSG_NOTICE([path recognized as directory])
#  else
#    AC_MSG_NOTICE([path not recognized as directory!: $cfsqp_path])
#  fi
  if test -f $cfsqp_path/cfsqp.c; then
    cp $cfsqp_path/cfsqp.c rigfilters/cfsqp/.
  fi
  if test -f $cfsqp_path/cfsqpusr.h; then
    cp $cfsqp_path/cfsqpusr.h rigfilters/cfsqp/.
  fi
  if test -f $cfsqp_path/qld.c; then
    cp $cfsqp_path/qld.c rigfilters/cfsqp/.
  fi
])
if test "$rigfilters" = 1; then
  if test -f rigfilters/cfsqp/cfsqp.c; then
    :
  else
    AC_MSG_NOTICE([CFSQP: Unable to locate cfsqp.c])
    rigfilters=0;
  fi
  if test -f rigfilters/cfsqp/cfsqpusr.h; then
    :
  else
    AC_MSG_NOTICE([CFSQP: Unable to locate cfsqpusr.h])
    rigfilters=0
  fi
  if test -f rigfilters/cfsqp/qld.c; then
    :
  else
    AC_MSG_NOTICE([CFSQP: Unable to locate qld.c])
    rigfilters=0
  fi
  if test "$rigfilters" = 0; then
    AC_MSG_NOTICE([CFSQP: Source not found; disabling rigorous filters])
    rigfilters=""
  else
    AC_MSG_NOTICE([CFSQP located, including rigfilters package in build])
    rigfilters="rigfilters"
  fi
fi
AC_SUBST(RIGFILTERS,$rigfilters)

# Default option for rigfilters is OFF
rigfilters=0;
# --with-rigfilters     install rigorous filters
AC_ARG_WITH(rigfilters,
[  --with-rigfilters       enable rigorous filters, requires CFSQP package],
[case $with_rigfilters in
   yes) AC_MSG_NOTICE([Rigorous filters selected])
	rigfilters=1
	;;
   no)  rigfilters=0
        ;;
    *)  echo "Ignoring unknown argument to --with-rigfilters: $with_rigfilters"
        ;;
esac])

#
# --with-mpi    MPI support for cmalign and cmsearch
AC_ARG_WITH(mpi,
[  --with-mpi            enable MPI for search and alignment],
[case $with_mpi in
   yes) AC_MSG_NOTICE([configured cmalign and cmsearch for MPI])
        AC_PROG_CC(mpicc)
        ;;
    no) ;;
     *) echo "Ignoring unknown argument to --with-mpi: $with_mpi"
        ;;	
  esac])


#
# end of optionally configured features
################################################################

#check for location of alloca() (used by cm2hmm)
AC_FUNC_ALLOCA

# Checks for system libraries
AC_CHECK_LIB(m, sqrt, LIBS="${LIBS} -lm")


# Configure external packages: the SQUID and Easel libraries.
#
AC_CONFIG_SUBDIRS(squid)
AC_CONFIG_SUBDIRS(easel)

# Write Infernal's config.h header
#
AC_CONFIG_HEADERS(src/config.h)

# Write out Makefiles in Infernal subdirectories.
AC_CONFIG_FILES(src/Makefile)
AC_CONFIG_FILES(testsuite/Makefile)
AC_CONFIG_FILES(rigfilters/Makefile)
AC_CONFIG_FILES(rigfilters/cfsqp/Makefile)
AC_CONFIG_FILES(rigfilters/cm2hmm-1.0/Makefile)


# remember that sneaky if statement that wraps the source
# configuration? well, it ends here.
else
  AC_MSG_NOTICE([This is a precompiled binary distribution - using abbreviated config])
fi

AC_CONFIG_FILES(Makefile)
AC_OUTPUT

