AM_CFLAGS = -I$(top_srcdir)/src
AM_CPPFLAGS = -I$(top_srcdir)/src

lib_LTLIBRARIES = librsmathutils.la
librsmathutils_la_SOURCES =                \
	correlation.c                          \
	filter.c                               \
	geom.c                                 \
	linalg.c                               \
	multivar_analysis.c                    \
	regression.c                           \
	stats.c                                \
	utils.c
librsmathutils_la_LIBADD = ../nifti/librsniftiutils.la @LIBOBJS@
