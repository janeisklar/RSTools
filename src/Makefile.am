if INCUBATOR
SUBDIRS = utils nifti maths . incubator batch
else
SUBDIRS = utils nifti maths . batch
endif

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\"
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\"

bin_PROGRAMS = rsbandpass rscorrelation rsinfo rsmask rsmotionscrubbing rsregression rsroi rstimecourse rsttest

rslibraries=utils/librsutils.a maths/librsmathutils.a nifti/librsniftiutils.a

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
rsroi_SOURCES = rsroi_ui.c rsroi_common.c rsroi.c
rstimecourse_SOURCES = rstimecourse_ui.c rstimecourse_common.c rstimecourse.c
rsttest_SOURCES = rsttest.c

lib_LIBRARIES = librstools.a
librstools_a_SOURCES = rsbandpass_ui.c rsbandpass_common.c rscorrelation_ui.c rscorrelation_common.c rsmotionscrubbing_ui.c rsmotionscrubbing_common.c rsregression_ui.c rsregression_common.c rsroi_ui.c rsroi_common.c rstimecourse_ui.c rstimecourse_common.c
librstools_a_LIBADD = @LIBOBJS@

pkginclude_HEADERS = rscommon.h maths/rsmathutils.h nifti/rsniftiutils.h utils/rsmemory.h utils/rsui.h utils/rsstring.h rsbandpass_ui.h rsbandpass_common.h rscorrelation_ui.h rscorrelation_common.h rsmotionscrubbing_ui.h rsmotionscrubbing_common.h rsregression_ui.h rsregression_common.h rsroi_ui.h rsroi_common.h rstimecourse_ui.h rstimecourse_common.h

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
