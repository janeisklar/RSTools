AUTOMAKE_OPTIONS = subdir-objects

lib_LTLIBRARIES = librsutils.la
librsutils_la_SOURCES = \
	rsmemory.c          \
	rsstring.c          \
	rsui.c
librsutils_la_LDFLAGS= -version-info 0:0:0
librsutils_la_LIBADD = @LIBOBJS@

pkginclude_HEADERS =    \
	rsmemory.h          \
	rsstring.h          \
	rsui.h