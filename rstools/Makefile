cc		=	gcc
CC		=	gcc
AR		=	ar
RANLIB  =	ranlib
CFLAGS	=	-ansi -pedantic -Woverlength-strings -g -std=gnu99 

FSLIO_INCS	=	-I../include
NIFTI_INCS	=	-I../include
ZNZ_INCS	=	-I/usr/include

GLIB_CXX            =   -static-libstdc++
FSLIO_LIBS	    =	-L../lib -lfslio
NIFTI_LIBS	    =	-L../lib -lniftiio
ZNZ_LIBS	    =	-L/usr/lib -L../lib -lznz -lm -lz
RSNIFTI_LIBS	=	-L. -lrsniftiutils
RSMATH_LIBS	    =	-L. -lrsmathutils

##############################################################
# R integration
R_HOME        = $(shell R RHOME)

ifeq ($(ARCH),X86_64) ## 64-bit Linux
	R_HOME        = $(shell /usr/local/rsTools/lib/R/bin/R RHOME)
endif

## include headers and libraries for R 
RCPPFLAGS     = $(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS      = $(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS         = $(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK       = $(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

## include headers and libraries for Rcpp interface classes
RCPPINCL      =	$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS      =	$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## include headers and libraries for RInside embedding classes
RINSIDEINCL   =	$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS   =	$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## compiler etc settings used in default make rules
CXX           = $(shell $(R_HOME)/bin/R CMD config CXX)
CPPFLAGS      =	-v -Wall $(shell $(R_HOME)/bin/R CMD config CPPFLAGS)
CXXFLAGS      = $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS)

##############################################################
# platform specific redefines (to use, set ARCH appropriately)

## ARCH = X86_64

DLL_EXT     = .dylib

ifeq ($(ARCH),SGI) ## SGI 32bit
ZNZ_INCS        =   -I/usr/freeware/include
ZNZ_LIBS        =   -L/usr/freeware/lib32 -L../lib -lznz -lm -lz
DLL_EXT         =   .so
else
ifeq ($(ARCH),I386) ## 32-bit Linux
ZNZ_INCS        =   -I/usr/include
ZNZ_LIBS        =   -L/usr/lib -L../lib -lznz -lm -lz
DLL_EXT         =   .so
GLIB_CXX        =   -L/usr/lib -lstdc++
else
ifeq ($(ARCH),X86_64) ## 64-bit Linux
ZNZ_INCS	=
ZNZ_LIBS	=   -L../lib -lznz -lm -lz
DLL_EXT         =   .so
GLIB_CXX        =   -static-libstdc++ ##-L/usr/lib -lstdc++
CXX             =   /usr/local/rsTools/lib/gcc/bin/g++-4.8
cc              =   /usr/local/rsTools/lib/gcc/bin/gcc-4.8
CC              =   /usr/local/rsTools/lib/gcc/bin/gcc-4.8
endif
endif
endif
LDLIBS        = $(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS) $(GLIB_CXX)

all:	mathutils niftiutils timecourse regression

clean:
	rm -f *.o *.a *.dylib

mathutils:
	@echo 
	@echo "***************************"
	@echo "*      RS Math Utils      *"
	@echo "***************************"
	$(CXX) $(CPPFLAGS) -c rsmathutils.cpp $(CXXFLAGS)
	$(CXX) $(CPPFLAGS) -shared rsmathutils.o -o librsmathutils$(DLL_EXT) $(CXXFLAGS) -shared-libgcc $(LDLIBS) ## -Wl,-Bstatic $(LDLIBS)

niftiutils:
	@echo
	@echo "***************************"
	@echo "*     RS Nifti Utils      *"
	@echo "***************************"
	$(CC) $(CFLAGS) -o rsniftiutils.o -c rsniftiutils.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS)
	$(AR) rcs librsniftiutils.a rsniftiutils.o
	$(RANLIB) librsniftiutils.a

timecourse:
	@echo
	@echo "***************************"
	@echo "*      RS Timecourse      *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rstimecourse.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS)
	$(CC) $(CFLAGS) rstimecourse.o -o rstimecourse $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS)

regression:
	@echo
	@echo "***************************"
	@echo "*      RS Regression      *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsregression.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS)
	##$(CC) $(CFLAGS) rsregression.o -o rsregression $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS)
	$(CXX) $(CPPFLAGS) rsregression.o -o rsregression $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(LDLIBS) $(RSMATH_LIBS)

help:
	@echo "all:      make the rstools"
	@echo "clean:    rm the rstools"
