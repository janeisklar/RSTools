AUTOMAKE_OPTIONS = subdir-objects
AM_CFLAGS = -DDATA_PATH=\"$(datadir)\" -DEXECUTABLES_PATH=\"$(bindir)\" -DCONFIG_PATH=\"$(sysconfdir)\" -I$(top_srcdir)/src
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\" -DEXECUTABLES_PATH=\"$(bindir)\" -DCONFIG_PATH=\"$(sysconfdir)\" -I$(top_srcdir)/src

lib_LTLIBRARIES = librsutils.la
librsutils_la_SOURCES = \
	rsio.c              \
	rsmemory.c          \
	rsstring.c          \
	rsui.c              \
	rsconfig.c
librsutils_la_LDFLAGS= -version-info 0:0:0
librsutils_la_LIBADD = @LIBOBJS@
