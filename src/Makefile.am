SUBDIRS = utils nifti maths . batch

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\"
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\"

rslibraries=utils/librsutils.la maths/librsmathutils.la nifti/librsniftiutils.la

LDADD=$(rslibraries)

bin_PROGRAMS =                                \
	rsbandpass                                \
	rscorrelation                             \
	rsinfo                                    \
	rsmask                                    \
	rsmotionscrubbing                         \
	rsregression                              \
	rsroi                                     \
	rstimecourse                              \
	rsttest
	
lib_LTLIBRARIES =                             \
	librsbandpass.la                          \
	librscorrelation.la                       \
	librsmotionscrubbing.la                   \
	librsregression.la                        \
	librsroi.la                               \
	librstimecourse.la

librsbandpass_la_SOURCES = rsbandpass_ui.c rsbandpass_common.c

librscorrelation_la_SOURCES = rscorrelation_ui.c rscorrelation_common.c
librsmotionscrubbing_la_SOURCES = rsmotionscrubbing_ui.c rsmotionscrubbing_common.c
librsregression_la_SOURCES = rsregression_ui.c rsregression_common.c 
librsroi_la_SOURCES = rsroi_ui.c rsroi_common.c
librstimecourse_la_SOURCES = rstimecourse_ui.c rstimecourse_common.c 

rsbandpass_SOURCES = rsbandpass.c
rsbandpass_LDADD = librsbandpass.la @LIBOBJS@ $(rslibraries)

rscorrelation_SOURCES = rscorrelation.c
rscorrelation_LDADD = librscorrelation.la @LIBOBJS@ $(rslibraries)

rsmotionscrubbing_SOURCES = rsmotionscrubbing.c
rsmotionscrubbing_LDADD = librsmotionscrubbing.la @LIBOBJS@ $(rslibraries)

rsregression_SOURCES = rsregression.c
rsregression_LDADD = librsregression.la @LIBOBJS@ $(rslibraries)

rsroi_SOURCES = rsroi.c
rsroi_LDADD = librsroi.la @LIBOBJS@ $(rslibraries)

rstimecourse_SOURCES = rstimecourse.c
rstimecourse_LDADD = librstimecourse.la @LIBOBJS@ $(rslibraries)

rsinfo_SOURCES = rsinfo.c
rsmask_SOURCES = rsmask.c
rsttest_SOURCES = rsttest.c

pkginclude_HEADERS =                       \
	rsbandpass_common.h                    \
	rsbandpass_ui.h                        \
	rscommon.h                             \
	rscorrelation_common.h                 \
	rscorrelation_ui.h                     \
	rsmotionscrubbing_common.h             \
	rsmotionscrubbing_ui.h                 \
	rsregression_common.h                  \
	rsregression_ui.h                      \
	rsroi_common.h                         \
	rsroi_ui.h                             \
	rstimecourse_common.h                  \
	rstimecourse_ui.h

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
