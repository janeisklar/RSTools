AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\" -I$(top_srcdir)/src
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\" -I$(top_srcdir)/src

rslibraries=                     \
        ../../utils/librsutils.la       \
        ../../maths/librsmathutils.la   \
        ../../nifti/librsniftiutils.la  \
        ../util/librsbatchutils.la

LDADD=$(rslibraries)

lib_LTLIBRARIES = \
	libRSBandpassPlugin.la \
	libRSRegressionPlugin.la \
	libRSCorrelationPlugin.la \
	libRSMotionScrubbingPlugin.la \
	libRSScrubbingPlugin.la \
	libRSRoiPlugin.la \
	libRSTimecoursePlugin.la \
	libUnixPlugin.la

libRSBandpassPlugin_la_SOURCES = rsbandpass/rsbandpass.cpp rsbandpass/tool/bandpass.cpp
libRSBandpassPlugin_la_LIBADD  = $(LDADD) ../../librsbandpass.la  @LIBOBJS@
libRSBandpassPlugin_la_LDFLAGS = -module -export-dynamic

libRSRegressionPlugin_la_SOURCES = rsregression/rsregression.cpp rsregression/tool/regression.cpp
libRSRegressionPlugin_la_LIBADD  = $(LDADD) ../../librsregression.la  @LIBOBJS@
libRSRegressionPlugin_la_LDFLAGS = -module -export-dynamic

libRSCorrelationPlugin_la_SOURCES = rscorrelation/rscorrelation.cpp rscorrelation/tool/correlation.cpp
libRSCorrelationPlugin_la_LIBADD  = $(LDADD) ../../librscorrelation.la  @LIBOBJS@
libRSCorrelationPlugin_la_LDFLAGS = -module -export-dynamic

libRSMotionScrubbingPlugin_la_SOURCES = rsmotionscrubbing/rsmotionscrubbing.cpp rsmotionscrubbing/tool/motionscrubbing.cpp
libRSMotionScrubbingPlugin_la_LIBADD  = $(LDADD) ../../librsmotionscrubbing.la  @LIBOBJS@
libRSMotionScrubbingPlugin_la_LDFLAGS = -module -export-dynamic

libRSScrubbingPlugin_la_SOURCES = rsscrubbing/rsscrubbing.cpp rsscrubbing/tool/scrubbing.cpp
libRSScrubbingPlugin_la_LIBADD  = $(LDADD) ../../librsscrubbing.la  @LIBOBJS@
libRSScrubbingPlugin_la_LDFLAGS = -module -export-dynamic

libRSRoiPlugin_la_SOURCES = rsroi/rsroi.cpp rsroi/tool/roi.cpp
libRSRoiPlugin_la_LIBADD  = $(LDADD) ../../librsroi.la  @LIBOBJS@
libRSRoiPlugin_la_LDFLAGS = -module -export-dynamic

libRSTimecoursePlugin_la_SOURCES = rstimecourse/rstimecourse.cpp rstimecourse/tool/timecourse.cpp
libRSTimecoursePlugin_la_LIBADD  = $(LDADD) ../../librstimecourse.la  @LIBOBJS@
libRSTimecoursePlugin_la_LDFLAGS = -module -export-dynamic

libUnixPlugin_la_SOURCES = unix/unix.cpp unix/tool/unix.cpp unix/task/unix.cpp
libUnixPlugin_la_LIBADD  = $(LDADD) ../util/librsbatchutils.la @LIBOBJS@
libUnixPlugin_la_LDFLAGS = -module -export-dynamic

pluginsdir = $(datadir)/$(PACKAGE)/plugins
plugins_DATA = \
	.libs/libRSBandpassPlugin.so \
	.libs/libRSRegressionPlugin.so \
	.libs/libRSCorrelationPlugin.so \
	.libs/libRSMotionScrubbingPlugin.so \
	.libs/libRSScrubbingPlugin.so \
	.libs/libRSRoiPlugin.so \
	.libs/libRSTimecoursePlugin.so \
	.libs/libUnixPlugin.so

jobpluginsdir = $(datadir)/$(PACKAGE)/jobs/plugins
jobplugins_DATA = \
	../../../share/jobs/plugins/rsbandpass.xsdext \
	../../../share/jobs/plugins/rsregression.xsdext \
	../../../share/jobs/plugins/rscorrelation.xsdext \
	../../../share/jobs/plugins/rsmotionscrubbing.xsdext \
	../../../share/jobs/plugins/rsscrubbing.xsdext \
	../../../share/jobs/plugins/rsroi.xsdext \
	../../../share/jobs/plugins/rstimecourse.xsdext \
	../../../share/jobs/plugins/unix.xsdext
