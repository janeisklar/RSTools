SUBDIRS = utils nifti maths . batch

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\" -DEXECUTABLES_PATH=\"$(bindir)\" -DCONFIG_PATH=\"$(sysconfdir)\" -I$(top_srcdir)/src
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\" -DEXECUTABLES_PATH=\"$(bindir)\" -DCONFIG_PATH=\"$(sysconfdir)\" -I$(top_srcdir)/src

rslibraries=utils/librsutils.la maths/librsmathutils.la nifti/librsniftiutils.la ../externals/librsfslio.la

LDADD=$(rslibraries)

bin_PROGRAMS =                                \
	rsbandpass                                \
	rscorrelation                             \
	rsdeoblique                               \
	rsinfo                                    \
	rsmask                                    \
	rsmotionscrubbing                         \
	rsscrubbing                               \
	rssmoothing                               \
	rsorientation                             \
	rsregression                              \
	rsresample                                \
	rsroi                                     \
	rstimecourse                              \
	rsttest                                   \
	rszeropadding                             \
	rsapplytransformation
	
lib_LTLIBRARIES =                             \
	librsbandpass.la                          \
	librscorrelation.la                       \
	librsdeoblique.la                         \
	librsinfo.la                              \
	librsmotionscrubbing.la                   \
	librsscrubbing.la                         \
	librssmoothing.la                         \
	librsregression.la                        \
	librsorientation.la                       \
	librsresample.la                          \
	librsroi.la                               \
	librstimecourse.la                        \
	librszeropadding.la                       \
	librsapplytransformation.la

librsbandpass_la_SOURCES = rsbandpass_ui.c rsbandpass_common.c
librsbandpass_la_LIBADD = @LIBOBJS@ $(rslibraries)
librscorrelation_la_SOURCES = rscorrelation_ui.c rscorrelation_common.c
librscorrelation_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsdeoblique_la_SOURCES = rsdeoblique_ui.c rsdeoblique_common.c
librsdeoblique_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsinfo_la_SOURCES = rsinfo_ui.c rsinfo_common.c
librsinfo_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsmotionscrubbing_la_SOURCES = rsmotionscrubbing_ui.c rsmotionscrubbing_common.c
librsmotionscrubbing_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsscrubbing_la_SOURCES = rsscrubbing_ui.c rsscrubbing_common.c
librsscrubbing_la_LIBADD = @LIBOBJS@ $(rslibraries)
librssmoothing_la_SOURCES = rssmoothing_ui.c rssmoothing_common.c
librssmoothing_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsorientation_la_SOURCES = rsorientation_ui.c rsorientation_common.c
librsorientation_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsregression_la_SOURCES = rsregression_ui.c rsregression_common.c 
librsregression_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsresample_la_SOURCES = rsresample_ui.c rsresample_common.c 
librsresample_la_LIBADD = @LIBOBJS@ $(rslibraries)
librsroi_la_SOURCES = rsroi_ui.c rsroi_common.c
librsroi_la_LIBADD = @LIBOBJS@ $(rslibraries)
librstimecourse_la_SOURCES = rstimecourse_ui.c rstimecourse_common.c
librstimecourse_la_LIBADD = @LIBOBJS@ $(rslibraries) 
librszeropadding_la_SOURCES = rszeropadding_ui.c rszeropadding_common.c
librszeropadding_la_LIBADD = @LIBOBJS@ $(rslibraries) 
librsapplytransformation_la_SOURCES = rsapplytransformation_ui.c rsapplytransformation_common.c
librsapplytransformation_la_LIBADD = librszeropadding.la @LIBOBJS@ $(rslibraries) 

rsbandpass_SOURCES = rsbandpass.c
rsbandpass_LDADD = librsbandpass.la @LIBOBJS@ $(rslibraries)

rscorrelation_SOURCES = rscorrelation.c
rscorrelation_LDADD = librscorrelation.la @LIBOBJS@ $(rslibraries)

rsdeoblique_SOURCES = rsdeoblique.c
rsdeoblique_LDADD = librsdeoblique.la @LIBOBJS@ $(rslibraries)

rsmotionscrubbing_SOURCES = rsmotionscrubbing.c
rsmotionscrubbing_LDADD = librsmotionscrubbing.la @LIBOBJS@ $(rslibraries)

rsscrubbing_SOURCES = rsscrubbing.c
rsscrubbing_LDADD = librsscrubbing.la @LIBOBJS@ $(rslibraries)

rssmoothing_SOURCES = rssmoothing.c
rssmoothing_LDADD = librssmoothing.la @LIBOBJS@ $(rslibraries)

rsorientation_SOURCES = rsorientation.c
rsorientation_LDADD = librsorientation.la @LIBOBJS@ $(rslibraries)

rsregression_SOURCES = rsregression.c
rsregression_LDADD = librsregression.la @LIBOBJS@ $(rslibraries)

rsresample_SOURCES = rsresample.c
rsresample_LDADD = librsresample.la @LIBOBJS@ $(rslibraries)

rsroi_SOURCES = rsroi.c
rsroi_LDADD = librsroi.la @LIBOBJS@ $(rslibraries)

rstimecourse_SOURCES = rstimecourse.c
rstimecourse_LDADD = librstimecourse.la @LIBOBJS@ $(rslibraries)

rsinfo_SOURCES = rsinfo.c
rsinfo_LDADD = librsinfo.la @LIBOBJS@ $(rslibraries)

rsmask_SOURCES = rsmask.c
rsttest_SOURCES = rsttest.c

rszeropadding_SOURCES = rszeropadding.c
rszeropadding_LDADD = librszeropadding.la @LIBOBJS@ $(rslibraries)

rsapplytransformation_SOURCES = rsapplytransformation.c
rsapplytransformation_LDADD = librsapplytransformation.la @LIBOBJS@ $(rslibraries)

nobase_pkginclude_HEADERS =                                   \
	rsbandpass_common.h                                       \
	rsbandpass_ui.h                                           \
	rscommon.h                                                \
	rscorrelation_common.h                                    \
	rscorrelation_ui.h                                        \
	rsdeoblique_common.h                                      \
	rsdeoblique_ui.h                                          \
	rsinfo_common.h                                           \
	rsinfo_ui.h                                               \
	rsmotionscrubbing_common.h                                \
	rsmotionscrubbing_ui.h                                    \
	rsorientation_common.h                                    \
	rsorientation_ui.h                                        \
	rsscrubbing_common.h                                      \
	rsscrubbing_ui.h                                          \
	rssmoothing_common.h                                      \
	rssmoothing_ui.h                                          \
	rsregression_common.h                                     \
	rsregression_ui.h                                         \
	rsresample_common.h                                       \
	rsresample_ui.h                                           \
	rsapplytransformation_common.h                            \
	rsapplytransformation_ui.h                                \
	rsroi_common.h                                            \
	rsroi_ui.h                                                \
	rstimecourse_common.h                                     \
	rstimecourse_ui.h                                         \
	rszeropadding_common.h                                    \
	rszeropadding_ui.h                                        \
	batch/plugins/rsbandpass/rsbandpass.hpp                   \
	batch/plugins/rsbandpass/tool/bandpass.hpp                \
	batch/plugins/rscorrelation/rscorrelation.hpp             \
	batch/plugins/rscorrelation/tool/correlation.hpp          \
	batch/plugins/rsmotionscrubbing/rsmotionscrubbing.hpp     \
	batch/plugins/rsmotionscrubbing/tool/motionscrubbing.hpp  \
	batch/plugins/rsscrubbing/rsscrubbing.hpp                 \
	batch/plugins/rsscrubbing/tool/scrubbing.hpp              \
	batch/plugins/rssmoothing/rssmoothing.hpp                 \
	batch/plugins/rssmoothing/tool/smoothing.hpp              \
	batch/plugins/rsorientation/rsorientation.hpp             \
	batch/plugins/rsorientation/tool/orientation.hpp          \
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
	batch/util/rsconfig.hpp                                   \
	batch/util/rsjob.hpp                                      \
	batch/util/rsjobparser.hpp                                \
	batch/util/rstask.hpp                                     \
	batch/util/rsunixtask.hpp                                 \
	batch/util/rsunixtool.hpp                                 \
	batch/util/rstool.hpp                                     \
	batch/util/rsunix.hpp                                     \
	maths/correlation.h                                       \
	maths/filter.h                                            \
	maths/geom.h                                              \
	maths/interpolation.h                                     \
	maths/linalg.h                                            \
	maths/multivar_analysis.h                                 \
	maths/regression.h                                        \
	maths/rsmathutils.h                                       \
	maths/stats.h                                             \
	maths/utils.h                                             \
	nifti/rsniftiutils.h                                      \
	nifti/headerinfo.h                                        \
	nifti/dicomfinder.h                                       \
	utils/rsio.h                                              \
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
	echo 'Libs: -L'${libdir}' -lrsbandpass -lrsbatchutils -lrscorrelation -lrsdeoblique -lrsfslio -lrsmathutils -lrsmotionscrubbing -lrsscrubbing -lrssmoothing -lrsorientation -lrsniftiutils -lrsregression -lrsresample -lrsroi -lrsinfo -lrstimecourse -lrszeropadding -lrsapplytransformation -lrsutils' >> $(srcdir)/rstools.pc
	echo 'Libs.private:' >> $(srcdir)/rstools.pc
	echo 'Cflags: -fopenmp -I'$(DESTDIR)$(includedir)'/rstools' >> $(srcdir)/rstools.pc
	
	$(MKDIR_P) $(DESTDIR)$(libdir)'/pkgconfig'
	$(INSTALL_PROGRAM_ENV) $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL_PROGRAM) $(srcdir)/rstools.pc $(DESTDIR)$(libdir)'/pkgconfig/rstools.pc' || exit $$?
