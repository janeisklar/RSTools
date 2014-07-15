include $(top_srcdir)/m4/autotroll.mk

SUBDIRS = util plugins .
AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = -DDATA_PATH=\"$(datadir)\"
AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\"

rslibraries=                     \
	../utils/librsutils.la       \
	../maths/librsmathutils.la   \
	../nifti/librsniftiutils.la  \
	util/librsbatchutils.la

bin_PROGRAMS = rsbatch rsjobeditor

LDADD=$(rslibraries)

rsbatch_SOURCES = rsbatch_common.cpp rsbatch_ui.cpp rsbatch.cpp

jobdir = $(datadir)/$(PACKAGE)/jobs
job_DATA = ../../share/jobs/job.xsd ../../share/jobs/dmn.job

rsjobeditor_SOURCES  = \
 jobeditor/ui/jobeditor.ui \
 rsjobeditor.cc \
 jobeditor/ui/ExtendedTabWidget.cpp                    jobeditor/ui/ExtendedTabWidget.moc.cpp \
 jobeditor/ui/ExtendedTabWidgetContainerExtension.cpp  jobeditor/ui/ExtendedTabWidgetContainerExtension.moc.cpp \
 jobeditor/ui/ExtendedTabWidgetExtensionFactory.cpp    jobeditor/ui/ExtendedTabWidgetExtensionFactory.moc.cpp \
 jobeditor/ui/ExtendedTabWidgetPlugin.cpp              jobeditor/ui/ExtendedTabWidgetPlugin.moc.cpp \
 jobeditor/rsjobeditorapplication.cpp                  jobeditor/rsjobeditorapplication.moc.cpp \
 jobeditor/rsjobeditorapplication.h
rsjobeditor_CXXFLAGS = $(QT_CXXFLAGS) $(AM_CXXFLAGS)
rsjobeditor_CPPFLAGS = $(QT_CPPFLAGS) $(AM_CPPFLAGS)
rsjobeditor_LDFLAGS  = $(QT_LDFLAGS) $(LDFLAGS)
rsjobeditor_LDADD    = $(QT_LIBS) $(LDADD)

BUILT_SOURCES = \
 jobeditor/ui/jobeditor.ui.h \
 jobeditor/ui/ExtendedTabWidget.moc.cpp \
 jobeditor/ui/ExtendedTabWidgetContainerExtension.moc.cpp \
 jobeditor/ui/ExtendedTabWidgetExtensionFactory.moc.cpp \
 jobeditor/ui/ExtendedTabWidgetPlugin.moc.cpp \
 jobeditor/rsjobeditorapplication.moc.cpp 
 
clean-local:
	find . -name '*.moc.cpp' -exec rm {} \;
	find . -name '*.ui.h' -exec rm {} \;