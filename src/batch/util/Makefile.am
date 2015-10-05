AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\" -DCONFIG_PATH=\"$(sysconfdir)\" -I$(top_srcdir)/src
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\" -DCONFIG_PATH=\"$(sysconfdir)\" -I$(top_srcdir)/src

rslibraries=../../utils/librsutils.a ../../librstools.a ../../maths/librsmathutils.a ../../nifti/librsniftiutils.a

lib_LTLIBRARIES = librsbatchutils.la
librsbatchutils_la_SOURCES =  \
 rsconfig.cpp \
 rsjob.cpp \
 rsjobparser.cpp \
 rstask.cpp \
 rstool.cpp \
 rsunix.cpp \
 rsunixtask.cpp \
 rsunixtool.cpp \
 plugin.cpp \
 pluginmanager.cpp
librsbatchutils_la_LDFLAGS=-shared
librsbatchutils_la_LIBADD = @LIBOBJS@ ../../utils/librsutils.la ../../nifti/librsniftiutils.la
