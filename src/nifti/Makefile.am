AM_CFLAGS = -I$(top_srcdir)/src
AM_CPPFLAGS = -I$(top_srcdir)/src

lib_LTLIBRARIES = librsniftiutils.la
librsniftiutils_la_SOURCES = rsniftiutils.c headerinfo.c dicomfinder.c
librsniftiutils_la_LIBADD = ../utils/librsutils.la ../../externals/librsfslio.la @LIBOBJS@
