AUTOMAKE_OPTIONS = foreign

# Allways put -I$(includedir) first to pick up includes from your own install area first
INCLUDES = -I$(includedir) -I$(ROOTSYS)/include -I$(OFFLINE_MAIN)/include

lib_LTLIBRARIES = \
  libhijing_analysis.la

noinst_HEADERS = \
  hijing_analysis.h

libhijing_analysis_la_SOURCES = \
  hijing_analysis.C \
  hijing_analysisDict.C

# Always put -L$(libdir) first to pick up libraries from your own install area first
# DLW: make sure that ROOTSYS appears before OFFLINE_MAIN to help pick up the root versions
# of the root libs instead of errant dummy OFFLINE versions.
libhijing_analysis_la_LIBADD = \
  -L$(libdir) \
  -L${ROOTSYS}/lib \
  -L$(OFFLINE_MAIN)/lib \
  -lfun4all

# Rule for generating CINT dictionaries from class headers.

hijing_analysisDict.C: \
  hijing_analysis.h    \
  hijing_analysisLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

