AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\"
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\"

lib_LTLIBRARIES = \
	libRSBandpassPlugin.la \
	libUnixPlugin.la
#libRSRegressionPlugin.la

libRSBandpassPlugin_la_SOURCES = rsbandpass/rsbandpass.cpp rsbandpass/tool/bandpass.cpp
libRSBandpassPlugin_la_LIBADD = ../util/librsbatchutils.la ../../librsbandpass.la  @LIBOBJS@
libRSBandpassPlugin_la_LDFLAGS = -module -export-dynamic

libUnixPlugin_la_SOURCES = unix/unix.cpp unix/tool/unix.cpp unix/task/unix.cpp
libUnixPlugin_la_LIBADD = ../util/librsbatchutils.la @LIBOBJS@
libUnixPlugin_la_LDFLAGS = -module -export-dynamic

#libRSRegressionPlugin_la_SOURCES = rsregression/rsregression.cpp
#libRSRegressionPlugin_la_LDFLAGS=-module -avoid-version -shared -export-dynamic
#libRSRegressionPlugin_la_LIBADD = @LIBOBJS@

pluginsdir = $(datadir)/$(PACKAGE)/plugins
plugins_DATA = \
	.libs/libRSBandpassPlugin.so \
	.libs/libUnixPlugin.so
#.libs/libRSRegressionPlugin.so

jobpluginsdir = $(datadir)/$(PACKAGE)/jobs/plugins
jobplugins_DATA = \
	../../../share/jobs/plugins/rsbandpass.xsdext \
	../../../share/jobs/plugins/unix.xsdext