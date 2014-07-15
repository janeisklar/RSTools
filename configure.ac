#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([RSTools], [0.9.0], [andre.hoffmann@meduniwien.ac.at])

ACLOCAL_AMFLAGS=-I m4
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([m4])

CFLAGS="-ansi -pedantic -Woverlength-strings -Wno-deprecated-register -g -std=c99 -fopenmp $CFLAGS"
CXXFLAGS="-ansi -pedantic -Woverlength-strings -Wno-deprecated-register -g -std=c++11 -fopenmp $CXXFLAGS"

# prepare for automake
AM_INIT_AUTOMAKE([foreign -Wall])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AM_PROG_AR
AC_PROG_RANLIB
#AC_LIBTOOL_DLOPEN
AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AC_PROG_INSTALL
LT_INIT

# Qt with AutoTroll.
AT_WITH_QT([+widgets +designer])

# determine git version information
RSTOOLS_VERSION_HASH=esyscmd([git --git-dir=./.git --work-tree=. rev-parse --short HEAD])
RSTOOLS_VERSION_DATE="esyscmd([git --git-dir=./.git --work-tree=. log -1 --pretty=format:'%ci'])"

AC_DEFINE_UNQUOTED(RSTOOLS_VERSION_HASH,
                   ["$RSTOOLS_VERSION_HASH"],
                   [A string containing the Git hash of the last commit]
)

AC_DEFINE_UNQUOTED(RSTOOLS_VERSION_DATE,
                   ["$RSTOOLS_VERSION_DATE"],
                   [A string containing the date of the last commit]
)

AC_DEFINE_UNQUOTED(GETTEXT_PACKAGE,
                   ["rstools"],
                   [The translation package for libgettext(unused)]
)

# Enabling/disabling the incubator
AC_ARG_ENABLE(
    [incubator],
    [  --enable-incubator      compile additional software that is not yet considered to be stable],
    [case "${enableval}" in
        yes) incubator=true ;;
        no)  incubator=false ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-incubator]) ;;
    esac],
    [incubator=false])
AM_CONDITIONAL([INCUBATOR], [test x$incubator = xtrue])

AM_COND_IF([INCUBATOR],
    [
        echo 'Incubator enabled'
        AC_CONFIG_FILES([src/incubator/Makefile])
    ],
    [
        echo 'Incubator disabled'
        AC_DEFINE(INCUBATOR_ENABLED, 0, [whether or not additional software that is not yet considered to be stable will be compiled])
    ]
)

# Check for FFTW3 library
AC_SEARCH_LIBS([fftw_execute], [fftw3],
    [AC_DEFINE( [RS_FFTW_ENABLED], [1], [Use FFTW])
       echo
       echo 'FFTW is being used'
       echo ] ,
    [AC_DEFINE( [RS_FFTW_ENABLED], [0], [Don't use FFTW])
       echo
       echo 'FFTW is not being used'
       echo ])

# Check for z-library
AC_SEARCH_LIBS([gzopen], [z],
    [] ,
    [  echo
       echo 'Lib z not found'
       echo 'Download from: http://www.zlib.net'
       echo 'or'
       echo 'apt-get install zlib1g'
       echo 
       exit 1 ])

# Check for m-library
AC_SEARCH_LIBS([sqrt], [m],
    [] ,
    [  echo
       echo 'C math library not found'
       echo 'apt-get install build-essential'
       echo 
       exit 1 ])
       
# Check for znz-library
AC_SEARCH_LIBS([znzopen], [znz],
    [] ,
    [  echo
       echo 'Lib znz not found'
       echo 'Download from: http://sourceforge.net/projects/niftilib/files/nifticlib/'
       echo 'or'
       echo 'apt-get install libnifti2 libnifti-dev'
       echo 
       exit 1 ])       

# Check for niftiio-library
AC_SEARCH_LIBS([nifti_image_load], [niftiio],
    [] ,
    [  echo
       echo 'Niftiio-library not found'
       echo 'Download from: http://sourceforge.net/projects/niftilib/files/nifticlib/'
       echo 'or'
       echo 'apt-get install libnifti2 libnifti-dev'
       echo 
       exit 1 ])

# Check for FSL library
AC_SEARCH_LIBS([FslGetDim], [fslio],
    [] ,
    [  echo
       echo 'Fslio not found'
       echo 'Download from: http://sourceforge.net/projects/niftilib/files/nifticlib/'
       echo 
       exit 1 ])
       
# Checking for Atlas
AC_ARG_ENABLE(
    [atlas],
    [  --disable-atlas         do not use ATLAS for BLAS-operations],
    [case "${enableval}" in
        yes) atlas=true ;;
        no)  atlas=false ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --disable-atlas]) ;;
    esac],
    [atlas=true])
AM_CONDITIONAL([ATLAS], [test x$atlas = xtrue])

AM_COND_IF([ATLAS],
    [
        AC_SEARCH_LIBS([cblas_ssymv], [atlas],
            [  echo 'Using ATLAS over GSL for BLAS-operations' 
               AC_DEFINE(USE_ATLAS, 1, [linked against ATLAS libraries])
            ] ,
            [  echo 'ATLAS not found, GSL will be used instead'
               echo 'Download from: http://math-atlas.sourceforge.net'
               echo 'or'
               echo 'apt-get install libatlas-base-dev'
               echo  
               AC_DEFINE(USE_ATLAS, 0, [linked against ATLAS libraries])
            ])
        AC_SEARCH_LIBS([ssymv], [cblas], [], [])
        AC_SEARCH_LIBS([pthread_create], [pthread], [], [])
    ],
    [
        echo 'Using GSL for BLAS-operations'
        AC_DEFINE(USE_ATLAS, 0, [linked against ATLAS libraries])
    ]
)

# Check for GSL
PKG_CHECK_MODULES([GSL], [gsl >= 1.15])
AC_SUBST([GSL_CFLAGS])
AC_SUBST([GSL_LIBS])
LIBS=["$GSL_LIBS $LIBS"]
CFLAGS=["$GSL_CFLAGS $CFLAGS"]
CXXFLAGS=["$GSL_CFLAGS $CXXFLAGS"]

# Check for GLIB
PKG_CHECK_MODULES([GLIB], [glib-2.0 >= 2.32.0])
AC_SUBST([GLIB_CFLAGS])
AC_SUBST([GLIB_LIBS])
LIBS=["$GLIB_LIBS $LIBS"]
CFLAGS=["$GLIB_CFLAGS $CFLAGS"]
CXXFLAGS=["$GLIB_CFLAGS $CXXFLAGS"]

# Check for GMODULE
PKG_CHECK_MODULES([GMODULE], [gmodule-2.0 >= 2.32.0])
AC_SUBST([GMODULE_CFLAGS])
AC_SUBST([GMODULE_LIBS])
LIBS=["$GMODULE_LIBS $LIBS"]
CFLAGS=["$GMODULE_CFLAGS $CFLAGS"]
CXXFLAGS=["$GMODULE_CFLAGS $CXXFLAGS"]

# Check for XERCES
PKG_CHECK_MODULES([XERCES], [xerces-c >= 3.0.0])
AC_SUBST([XERCES_CFLAGS])
AC_SUBST([XERCES_LIBS])
LIBS=["$XERCES_LIBS $LIBS"]
CXXFLAGS=["$XERCES_CFLAGS $CXXFLAGS"]

# Checks for header files.
AC_CHECK_HEADERS([float.h string.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

# Define included makefiles
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
		 src/Makefile
 		 src/maths/Makefile
 		 src/nifti/Makefile
 		 src/batch/Makefile
		 src/batch/plugins/Makefile
		 src/batch/util/Makefile
 		 src/utils/Makefile])

AC_OUTPUT
