SUBDIRS = nifti maths .

bin_PROGRAMS = rsbandpass rscentrality rscorrelation rsfastcentrality rsfillholes rsinfo rsmask rsmotionscrubbing rsregression rsroi rstimecourse rsttest

rslibraries=librsmathutils.a librsniftiutils.a

LDADD=$(rslibraries)

rsbandpass_SOURCES = rsbandpass_common.c rsbandpass.c
#rsbandpass_LDADD = $(LIBOBJS)
#rsbandpass_LDADD = -Lnifti -Lmaths -lrsmathutils -lrsniftiutils
rscentrality_SOURCES = rscentrality.c
rscorrelation_SOURCES = rscorrelation_common.c rscorrelation.c
rsfastcentrality_SOURCES = rsfastcentrality.c
rsfillholes_SOURCES = rsfillholes.c
rsinfo_SOURCES = rsinfo.c
rsmask_SOURCES = rsmask.c
rsmotionscrubbing_SOURCES = rsmotionscrubbing.c
rsmotionscrubbing_LDADD = rsregression_common.o @LIBOBJS@ $(rslibraries)
rsregression_SOURCES = rsregression.c rsregression_common.c
rsroi_SOURCES = rsroi.c
rstimecourse_SOURCES = rstimecourse.c
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
