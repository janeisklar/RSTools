SUBDIRS = utils nifti maths . batch

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\" -I$(top_srcdir)/src
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\" -I$(top_srcdir)/src

rslibraries=utils/librsutils.la maths/librsmathutils.la nifti/librsniftiutils.la ../externals/librsfslio.la

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
librsbandpass_la_LIBADD = @LIBOBJS@ $(rslibraries)
librscorrelation_la_SOURCES = rscorrelation_ui.c rscorrelation_common.c
librscorrelation_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsmotionscrubbing_la_SOURCES = rsmotionscrubbing_ui.c rsmotionscrubbing_common.c
librsmotionscrubbing_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsregression_la_SOURCES = rsregression_ui.c rsregression_common.c 
librsregression_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsroi_la_SOURCES = rsroi_ui.c rsroi_common.c
librsroi_la_LIBADD = @LIBOBJS@ $(rslibraries)
librstimecourse_la_SOURCES = rstimecourse_ui.c rstimecourse_common.c
librstimecourse_la_LIBADD = @LIBOBJS@ $(rslibraries) 

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

nobase_pkginclude_HEADERS =                                   \
	rsbandpass_common.h                                       \
	rsbandpass_ui.h                                           \
	rscommon.h                                                \
	rscorrelation_common.h                                    \
	rscorrelation_ui.h                                        \
	rsmotionscrubbing_common.h                                \
	rsmotionscrubbing_ui.h                                    \
	rsregression_common.h                                     \
	rsregression_ui.h                                         \
	rsroi_common.h                                            \
	rsroi_ui.h                                                \
	rstimecourse_common.h                                     \
	rstimecourse_ui.h                                         \
	batch/plugins/rsbandpass/rsbandpass.hpp                   \
	batch/plugins/rsbandpass/tool/bandpass.hpp                \
	batch/plugins/rscorrelation/rscorrelation.hpp             \
	batch/plugins/rscorrelation/tool/correlation.hpp          \
	batch/plugins/rsmotionscrubbing/rsmotionscrubbing.hpp     \
	batch/plugins/rsmotionscrubbing/tool/motionscrubbing.hpp  \
	batch/plugins/rsregression/rsregression.hpp               \
	batch/plugins/rsregression/tool/regression.hpp            \
	batch/plugins/rsroi/rsroi.hpp                             \
	batch/plugins/rsroi/tool/roi.hpp                          \
	batch/plugins/rstimecourse/rstimecourse.hpp               \
	batch/plugins/rstimecourse/tool/timecourse.hpp            \
	batch/plugins/unix/task/unix.hpp                          \
	batch/plugins/unix/tool/unix.hpp                          \
	batch/plugins/unix/unix.hpp                               \
	batch/rsbatch_common.hpp                                  \
	batch/rsbatch_ui.hpp                                      \
	batch/util/outputCatcher.hpp                              \
	batch/util/parseErrorHandler.hpp                          \
	batch/util/plugin.hpp                                     \
	batch/util/pluginmanager.hpp                              \
	batch/util/rsjob.hpp                                      \
	batch/util/rsjobparser.hpp                                \
	batch/util/rstask.hpp                                     \
	batch/util/rstool.hpp                                     \
	maths/correlation.h                                       \
	maths/filter.h                                            \
	maths/geom.h                                              \
	maths/linalg.h                                            \
	maths/multivar_analysis.h                                 \
	maths/regression.h                                        \
	maths/rsmathutils.h                                       \
	maths/stats.h                                             \
	maths/utils.h                                             \
	nifti/rsniftiutils.h                                      \
	utils/rsmemory.h                                          \
	utils/rsstring.h                                          \
	utils/rsui.h                                              

install-exec-hook:			   
	echo 'prefix='$(prefix) > $(srcdir)/rstools.pc
	echo 'exec_prefix='$(DESTDIR)$(bindir) >> $(srcdir)/rstools.pc
	echo 'libdir='$(DESTDIR)$(libdir) >> $(srcdir)/rstools.pc
	echo 'includedir='$(DESTDIR)$(includedir) >> $(srcdir)/rstools.pc
	echo '' >> $(srcdir)/rstools.pc
	echo 'Name: RSTools' >> $(srcdir)/rstools.pc
	echo 'Description: Resting-state toolbox' >> $(srcdir)/rstools.pc
	echo 'Requires: glib-2.0 >= 2.32.0, gsl >= 1.15, gmodule-2.0 >= 2.32.0, xerces-c >= 3.0.0' >> $(srcdir)/rstools.pc
	echo 'Version: '$(VERSION) >> $(srcdir)/rstools.pc
	echo 'Libs: -L'${libdir}' -lrsbandpass -lrsbatchutils -lrscorrelation -lrsfslio -lrsmathutils -lrsmotionscrubbing -lrsniftiutils -lrsregression -lrsroi -lrstimecourse -lrsutils' >> $(srcdir)/rstools.pc
	echo 'Libs.private:' >> $(srcdir)/rstools.pc
	echo 'Cflags: -fopenmp -I'$(DESTDIR)$(includedir)'/rstools' >> $(srcdir)/rstools.pc
	
	$(MKDIR_P) $(DESTDIR)$(libdir)'/pkgconfig'
	$(INSTALL_PROGRAM_ENV) $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL_PROGRAM) $(srcdir)/rstools.pc $(DESTDIR)$(libdir)'/pkgconfig/rstools.pc' || exit $$?
