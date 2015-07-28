AUTOMAKE_OPTIONS = foreign

CINTINCLUDES = -pthread -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I$(ROOTSYS)/root/include

INCLUDES = $(CINTINCLUDES) 

lib_LTLIBRARIES = \
  libCombinedSimple.la 

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
nobase_dist_mydata_DATA = \
 inputs/run10/AA_trig_eff_0_v0.root \
 inputs/run10/AA_trig_eff_1_v0.root \
 inputs/run10/AA_trig_eff_2_v0.root \
 inputs/run10/AA_trig_eff_3_v0.root \
 inputs/run10/sharkfin_projection.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent0_fg_sys0.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent0_fg_sys1.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent0_fg_sys-1.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent0_bg_sys0.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent0_bg_sys1.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent0_bg_sys-1.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent1_fg_sys0.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent1_fg_sys1.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent1_fg_sys-1.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent1_bg_sys0.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent1_bg_sys1.root \
 inputs/run10/filltimecorrs_PisaEffcorr_cent1_bg_sys-1.root \
 inputs/run10/badtowers.dat \
 inputs/run11/badtowers.dat

libCombinedSimple_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib\
  -lfun4all \
  -lSubsysReco \
  -lnanoDST \
  -llvl2 \
  -lfun4allfuncs \
  -lert \
  -lrecal \
  -lTHmul \
  -lFROG \
  -lTOAD

libCombinedSimple_la_SOURCES = \
  CombinedSimple.C \
  combOnetrigbgSimple.C \
  combOnetrigbgTrack.C \
  CombinedSimple_Dict.C

 BUILT_SOURCES = \
   testexternals.C

 noinst_PROGRAMS = \
   testexternals

 testexternals_LDADD = \
   libCombinedSimple.la

 testexternals.C:
	 echo "//*** this is a generated file. Do not commit, do not edit" > $@
	 echo "int main()" >> $@
	 echo "{" >> $@
	 echo "  return 0;" >> $@
	 echo "}" >> $@


CombinedSimple_Dict.C: CombinedSimple.h combOnetrigbgSimple.h combOnetrigbgTrack.h CombinedSimpleLinkDef.h 
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(CINTINCLUDES) $^

clean-local:
	rm -f *Dict*

