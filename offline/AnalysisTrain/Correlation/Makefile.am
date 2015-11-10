AUTOMAKE_OPTIONS = foreign

ROOTCFLAGS = $(shell root-config --cflags)

INCLUDES = -I$(includedir) \
-I$(OFFLINE_MAIN)/include \
-I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libCorrelation.la

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
nobase_dist_mydata_DATA = \
  Warnmap/d_Au_warnmap.txt \
  Warnmap/d_Au_warnmap_0to5GeV_cut.txt \
  Warnmap/d_Au_warnmap_above5GeV_cut.txt \
  Warnmap/newwarnmap_0to5GeV.txt \
  Warnmap/newwarnmap_above5GeV.txt \
  Warnmap/extrawarnmap.txt \
  Warnmap/warnmap_0to5_run10.txt \
  Warnmap/warnmap_above5_run10.txt \
  Warnmap/extrawarnmap_run10.txt \
  inputs/run10/AA_trig_eff_0_v0.root \
  inputs/run10/AA_trig_eff_1_v0.root \
  inputs/run10/AA_trig_eff_2_v0.root \
  inputs/run10/AA_trig_eff_3_v0.root \
  inputs/run10/sharkfin_projection.root \
  inputs/run10/filltimecorrs_PisaEffcorr_cent0_bg_sys0.root \
  inputs/run11/AA_11_trig_eff_0_v0.root \
  inputs/run11/AA_11_trig_eff_1_v0.root \
  inputs/run11/AA_11_trig_eff_2_v0.root \
  inputs/run11/AA_11_trig_eff_3_v0.root \
  inputs/run11/filltimeweighting.root \
  inputs/run11/run11_hadron_eff_func.root \
  inputs/run11/v2_inputs.root \
  inputs/run8/sharkfin_miss.root \
  inputs/run8/dA_run8_trig_eff_0_v0.root \
  inputs/run8/dA_run8_trig_eff_1_v0.root \
  inputs/run8/dA_run8_trig_eff_2_v0.root \
  inputs/run8/dA_run8_trig_eff_3_v0.root \
  inputs/run8/chhadron_eff_dAu_C0.root \
  inputs/run8/filltimeweighting.root

libCorrelation_la_LIBADD = \
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

noinst_HEADERS = \
  CorrelationLinkDef.h

include_HEADERS = \
  AParticle.h \
  ACluster.h \
  ATrack.h \
  APiZero.h \
  AMixingTree.h \
  CorrelationFunctions.h
  
libCorrelation_la_SOURCES = \
  AParticle.C \
  ACluster.C \
  ATrack.C \
  APiZero.C \
  AEvent.C \
  AMixingPool.C \
  AMixingTree.C \
  Warnmap.C \
  Correlation.C

nodist_libCorrelation_la_SOURCES = Correlation_Dict.C

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_SOURCES = testexternals.C

testexternals_LDADD = \
  libCorrelation.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

Correlation_Dict.C: AParticle.h ACluster.h ATrack.h APiZero.h AEvent.h AMixingPool.h AMixingTree.h Warnmap.h Correlation.h CorrelationLinkDef.h
	rootcint -f $@ -c -O -Wall -fPIC $(ROOTCFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f ${BUILT_SOURCES}
	rm -f *Dict* G__*LinkDef.h
