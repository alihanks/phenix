AUTOMAKE_OPTIONS = foreign

ROOTCFLAGS = $(shell root-config --cflags)

INCLUDES = -I$(includedir) \
-I$(OFFLINE_MAIN)/include \
-I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libCorrelationPlots.la

libCorrelationPlots_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lTOAD \
  @ROOTLIBS@ \
  -lPhysics \
  -lfun4all \
  -lSubsysReco \
  -lrecal \
  -lfun4allfuncs \
  -lTHmul

include_HEADERS = \
  MakeCFs.h \
  MakeJFs.h \
  MakeCombinedHistos.h \
  MakeDir.h \
  MakeWeightedJFs.h

noinst_HEADERS = MakePlotsLinkDef.h

libCorrelationPlots_la_SOURCES = \
  MakeCFs.C \
  MakeJFs.C \
  MakeCombinedHistos.C \
  MakeDir.C \
  MakeWeightedJFs.C

nodist_libCorrelationPlots_la_SOURCES = CorrelationPlots_Dict.C

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_SOURCES = testexternals.C

testexternals_LDADD = \
  libCorrelationPlots.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

CorrelationPlots_Dict.C: MakeCFs.h MakeJFs.h MakeDir.h MakeWeightedJFs.h MakeCombinedHistos.h MakePlotsLinkDef.h
	rootcint -f $@ -c -O -Wall -fPIC $(ROOTCFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f ${BUILT_SOURCES}
	rm -f *Dict* G__*LinkDef.h
