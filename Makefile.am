ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = -DDATA_PATH=\"$(datadir)\"
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\"
AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip
EXTRA_DIST = autogen.sh

SUBDIRS = externals src

nobase_pkginclude_HEADERS =      \
	externals/fslio/fslio.h      \
	externals/fslio/dbh.h

sysconf_DATA=etc/rstools.conf