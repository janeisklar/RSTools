if INCUBATOR
SUBDIRS = nifti maths . incubator
else
SUBDIRS = nifti maths .
endif

bin_PROGRAMS = rsbandpass rscorrelation rsinfo rsmask rsmotionscrubbing rsregression rsroi rstimecourse rsttest

rslibraries=maths/librsmathutils.a nifti/librsniftiutils.a

LDADD=$(rslibraries)

rsbandpass_SOURCES = rsbandpass_ui.c rsbandpass_common.c rsbandpass.c
#rsbandpass_LDADD = $(LIBOBJS)
#rsbandpass_LDADD = -Lnifti -Lmaths -lrsmathutils -lrsniftiutils
rscorrelation_SOURCES = rscorrelation_ui.c rscorrelation_common.c rscorrelation.c 
rsinfo_SOURCES = rsinfo.c
rsmask_SOURCES = rsmask.c
rsmotionscrubbing_SOURCES = rsmotionscrubbing_ui.c rsmotionscrubbing_common.c rsmotionscrubbing.c
#rsmotionscrubbing_LDADD = rsregression_common.o @LIBOBJS@ $(rslibraries)
rsregression_SOURCES = rsregression_ui.c rsregression_common.c rsregression.c
rsroi_SOURCES = rsroi.c
rstimecourse_SOURCES = rstimecourse_ui.c rstimecourse_common.c rstimecourse.c
rsttest_SOURCES = rsttest.c

pkginclude_HEADERS = rscommon.h maths/rsmathutils.h nifti/rsniftiutils.h

install-exec-hook:
	ln -fs $(DESTDIR)$(bindir)/rsbandpass$(EXEEXT) \
	       $(DESTDIR)$(bindir)/rsbandpass2$(EXEEXT)
	ln -fs $(DESTDIR)$(bindir)/rsregression$(EXEEXT) \
	       $(DESTDIR)$(bindir)/rsregression2$(EXEEXT)
	ln -fs $(DESTDIR)$(bindir)/rsmotionscrubbing$(EXEEXT) \
	       $(DESTDIR)$(bindir)/rsmotionscrubbing2$(EXEEXT)
	ln -fs $(DESTDIR)$(bindir)/rsttest$(EXEEXT) \
	       $(DESTDIR)$(bindir)/rsttest2$(EXEEXT)
	ln -fs $(DESTDIR)$(bindir)/rstimecourse$(EXEEXT) \
	       $(DESTDIR)$(bindir)/rstimecourse2$(EXEEXT)
	ln -fs $(DESTDIR)$(bindir)/rscorrelation$(EXEEXT) \
               $(DESTDIR)$(bindir)/rscorrelation2$(EXEEXT)
