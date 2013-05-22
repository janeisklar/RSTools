CC	=	gcc
CFLAGS	=	-ansi -pedantic -Woverlength-strings -g

FSLIO_INCS	=	-I../include
NIFTI_INCS	=	-I../include
ZNZ_INCS	=	-I/usr/include

FSLIO_LIBS	=	-L../lib -lfslio
NIFTI_LIBS	=	-L../lib -lniftiio
ZNZ_LIBS	=	-L/usr/lib -L../lib -lznz -lm -lz


##############################################################
# platform specific redefines (to use, set ARCH appropriately)

## ARCH = X86_64

ifeq ($(ARCH),SGI) ## SGI 32bit
ZNZ_INCS	=	-I/usr/freeware/include
ZNZ_LIBS	=	-L/usr/freeware/lib32 -L../lib -lznz -lm -lz
else
ifeq ($(ARCH),I386) ## 32-bit Linux
ZNZ_INCS	=	-I/usr/include
ZNZ_LIBS	=	-L/usr/lib -L../lib -lznz -lm -lz
else
ifeq ($(ARCH),X86_64) ## 64-bit Linux
ZNZ_INCS	=
ZNZ_LIBS	=	-L../lib -lznz -lm -lz
endif
endif
endif



all:	timecourse

clean:
	rm -f *.o *.a

timecourse: 
	$(CC) $(CFLAGS) -o rstimecourse rstimecourse.c $(FSLIO_INCS) $(NIFTI_INCS) $(ZNZ_INCS) $(FSLIO_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS)

help:
	@echo "all:      make the rstools"
	@echo "clean:    rm the rstools"
