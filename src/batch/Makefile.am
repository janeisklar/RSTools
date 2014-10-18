SUBDIRS = util plugins .
AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\" -I$(top_srcdir)/src
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\" -I$(top_srcdir)/src

rslibraries=                     \
	../utils/librsutils.la       \
	../maths/librsmathutils.la   \
	../nifti/librsniftiutils.la  \
	util/librsbatchutils.la

bin_PROGRAMS = rsbatch

LDADD=$(rslibraries)

rsbatch_SOURCES = rsbatch_common.cpp rsbatch_ui.cpp rsbatch.cpp

jobdir = $(datadir)/$(PACKAGE)/jobs
job_DATA = ../../share/jobs/job.xsd ../../share/jobs/dmn.job ../../share/jobs/empty.job
