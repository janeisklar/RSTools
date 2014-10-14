AUTOMAKE_OPTIONS = subdir-objects
AM_CFLAGS = -I$(top_srcdir)/src
AM_CPPFLAGS = -I$(top_srcdir)/src

lib_LTLIBRARIES = librsutils.la
librsutils_la_SOURCES = \
	rsio.c              \
	rsmemory.c          \
	rsstring.c          \
	rsui.c
librsutils_la_LDFLAGS= -version-info 0:0:0
librsutils_la_LIBADD = @LIBOBJS@
