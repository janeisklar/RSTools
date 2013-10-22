cc		=	gcc
CC		=	gcc
AR		=	ar
RANLIB  =	ranlib

FSLIO_INCS	=	-I../include
NIFTI_INCS	=	-I../include
ZNZ_INCS	=	-I/usr/include

GLIB_CXX            =   -static-libstdc++
FSLIO_LIBS	    =	-L../lib -lfslio
NIFTI_LIBS	    =	-L../lib -lniftiio
ZNZ_LIBS	    =	-L/usr/lib -L../lib -lznz -lm -lz
RSNIFTI_LIBS	=	-L. -lrsniftiutils

ifndef NOFFTW
FFTW_LIBS       =   -lfftw3
FFTW_CFLAGS     = 	-DRS_FFTW_ENABLED=1
else
FFTW_LIBS	    =
FFTW_CFLAGS     = 	
endif

ifdef NATIVE
NATIVE_CFLAGS   =   -march=native
endif

CFLAGS          =   -ansi -pedantic -Woverlength-strings -g -std=gnu99 $(FFTW_CFLAGS) $(NATIVE_CFLAGS)
RSMATH_LIBS	    =	-L. -lrsmathutils -lgsl -lgslcblas $(FFTW_LIBS)

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
#AR              =   /usr/local/rsTools/lib/gcc/bin/gcc-ar-4.8
#RANLIB          =   /usr/local/rsTools/lib/gcc/bin/gcc-ranlib-4.8
endif
endif
endif
LDLIBS        = $(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS) $(GLIB_CXX) -lgsl -lgslcblas

all:	mathutils niftiutils timecourse regression bandpass roi correlation centrality fillholes fastcentrality ttest mask motionscrubbing

clean:
	rm -f *.o *.a *.dylib ../lib/librs*
	rm -f rsregression
	rm -f rsregression2
	rm -f rstimecourse
	rm -f rstimecourse2
	rm -f rsroi
	rm -f rscorrelation
	rm -f rscorrelation2
	rm -f rsbandpass
	rm -f rsbandpass2
	rm -f rscentrality2
	rm -f rsfillholes2
	rm -f rsfastcentrality2
	rm -f rsttest2
	rm -f rsmask2
	rm -f rsmotionscrubbing2

mathutils:
	@echo 
	@echo "***************************"
	@echo "*      RS Math Utils      *"
	@echo "***************************"
	$(CC) $(CFLAGS) -o rsmathutils.o -c rsmathutils.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(AR) rcs librsmathutils.a rsmathutils.o
	$(RANLIB) librsmathutils.a

niftiutils:
	@echo
	@echo "***************************"
	@echo "*     RS Nifti Utils      *"
	@echo "***************************"
	$(CC) $(CFLAGS) -o rsniftiutils.o -c rsniftiutils.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(AR) rcs librsniftiutils.a rsniftiutils.o
	$(RANLIB) librsniftiutils.a

timecourse:
	@echo
	@echo "***************************"
	@echo "*      RS Timecourse      *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rstimecourse.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rstimecourse.o -o rstimecourse $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) -fopenmp
	$(CC) $(CFLAGS) -c rstimecourse2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rstimecourse2.o -o rstimecourse2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) -fopenmp

regression:
	@echo
	@echo "***************************"
	@echo "*      RS Regression      *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsregression_common.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) -c rsregression.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) -c rsregression2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rsregression_common.o rsregression.o -o rsregression $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp
	$(CC) $(CFLAGS) rsregression_common.o rsregression2.o -o rsregression2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp

bandpass:
	@echo
	@echo "***************************"
	@echo "*       RS Bandpass       *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsbandpass_common.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) -c rsbandpass.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) -c rsbandpass2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rsbandpass_common.o rsbandpass.o -o rsbandpass $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp
	$(CC) $(CFLAGS) rsbandpass_common.o rsbandpass2.o -o rsbandpass2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp

roi:
	@echo
	@echo "***************************"
	@echo "*         RS ROI          *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsroi.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rsroi.o -o rsroi $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp

correlation:
	@echo
	@echo "***************************"
	@echo "*     RS Correlation      *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rscorrelation.c -c rscorrelation_common.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) -c rscorrelation2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rscorrelation_common.o rscorrelation.o -o rscorrelation $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp
	$(CC) $(CFLAGS) rscorrelation_common.o rscorrelation2.o -o rscorrelation2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp

centrality:
	@echo
	@echo "***************************"
	@echo "*     RS Centrality       *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rscentrality2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rscentrality2.o -o rscentrality2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp

fillholes:
	@echo
	@echo "***************************"
	@echo "*     RS Fill holes       *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsfillholes2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rsfillholes2.o -o rsfillholes2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp

fastcentrality:
	@echo
	@echo "***************************"
	@echo "*   RS Fast Centrality    *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsfastcentrality2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rsfastcentrality2.o -o rsfastcentrality2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp
    
ttest:
	@echo
	@echo "***************************"
	@echo "*        RS T-Test        *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsttest2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rsttest2.o -o rsttest2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp
    
mask:
	@echo
	@echo "***************************"
	@echo "*        RS Mask          *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsmask2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rsmask2.o -o rsmask2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp

motionscrubbing:
	@echo
	@echo "***************************"
	@echo "*  RS Motion Scrbubbing   *"
	@echo "***************************"
	$(CC) $(CFLAGS) -c rsmotionscrubbing2.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) -fopenmp
	$(CC) $(CFLAGS) rsregression_common.o rsmotionscrubbing2.o -o rsmotionscrubbing2 $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS) $(RSNIFTI_LIBS) $(RSMATH_LIBS) -fopenmp

help:
	@echo "all:      make the rstools"
	@echo "clean:    rm the rstools"
