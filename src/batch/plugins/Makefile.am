AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\"
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\"

lib_LTLIBRARIES = \
	libRSBandpassPlugin.la \
	libRSRegressionPlugin.la \
	libRSCorrelationPlugin.la \
	libUnixPlugin.la

libRSBandpassPlugin_la_SOURCES = rsbandpass/rsbandpass.cpp rsbandpass/tool/bandpass.cpp
libRSBandpassPlugin_la_LIBADD  = ../util/librsbatchutils.la ../../librsbandpass.la  @LIBOBJS@
libRSBandpassPlugin_la_LDFLAGS = -module -export-dynamic

libRSRegressionPlugin_la_SOURCES = rsregression/rsregression.cpp rsregression/tool/regression.cpp
libRSRegressionPlugin_la_LIBADD  = ../util/librsbatchutils.la ../../librsregression.la  @LIBOBJS@
libRSRegressionPlugin_la_LDFLAGS = -module -export-dynamic

libRSCorrelationPlugin_la_SOURCES = rscorrelation/rscorrelation.cpp rscorrelation/tool/correlation.cpp
libRSCorrelationPlugin_la_LIBADD  = ../util/librsbatchutils.la ../../librscorrelation.la  @LIBOBJS@
libRSCorrelationPlugin_la_LDFLAGS = -module -export-dynamic

libUnixPlugin_la_SOURCES = unix/unix.cpp unix/tool/unix.cpp unix/task/unix.cpp
libUnixPlugin_la_LIBADD  = ../util/librsbatchutils.la @LIBOBJS@
libUnixPlugin_la_LDFLAGS = -module -export-dynamic

pluginsdir = $(datadir)/$(PACKAGE)/plugins
plugins_DATA = \
	.libs/libRSBandpassPlugin.so \
	.libs/libRSRegressionPlugin.so \
	.libs/libRSCorrelationPlugin.so \
	.libs/libUnixPlugin.so

jobpluginsdir = $(datadir)/$(PACKAGE)/jobs/plugins
jobplugins_DATA = \
	../../../share/jobs/plugins/rsbandpass.xsdext \
	../../../share/jobs/plugins/rsregression.xsdext \
	../../../share/jobs/plugins/rscorrelation.xsdext \
	../../../share/jobs/plugins/unix.xsdext