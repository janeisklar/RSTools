#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([RSTools], [0.8.1], [andre.hoffmann@meduniwien.ac.at])

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

CFLAGS="-ansi -pedantic -Woverlength-strings -g -std=gnu99 -fopenmp $CFLAGS"

# prepare for automake
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/rsbandpass.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
		 src/Makefile
 		 src/maths/Makefile
 		 src/nifti/Makefile])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_LIBTOOL
AC_PROG_INSTALL

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
    [  --disable-atlas    Do not use ATLAS for BLAS-operations],
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

AC_OUTPUT
