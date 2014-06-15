AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\"
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\"

lib_LIBRARIES = librsjob.a
librsjob_a_SOURCES = rsjob.cpp rsjobparser.cpp rstask.cpp execution/timecourse.cpp execution/regression.cpp execution/bandpass.cpp execution/motionscrubbing.cpp execution/correlation.cpp execution/roi.cpp execution/unix.cpp
librsjob_a_LIBADD = @LIBOBJS@

bin_PROGRAMS = rsbatch

rslibraries=../maths/librsmathutils.a ../nifti/librsniftiutils.a ../librstools.a
LDADD=$(rslibraries) librsjob.a

rsbatch_SOURCES = rsbatch_common.cpp rsbatch_ui.cpp rsbatch.cpp

jobdir = $(datadir)/$(PACKAGE)/jobs
job_DATA = ../../share/jobs/job.xsd ../../share/jobs/dmn.job