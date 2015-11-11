# void draw_isoDist_EPdep_sumIso(TString outName="AllQCDPhoton30", bool isPromptPho=1, TString treePath="ggHiNtuplizer",int rad=2)

# AllQCDPhoton30 , ggHiNtuplizer
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("AllQCDPhoton30", 1, "ggHiNtuplizer", 2)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("AllQCDPhoton30", 1, "ggHiNtuplizer", 3)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("AllQCDPhoton30", 1, "ggHiNtuplizer", 4)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("AllQCDPhoton30", 1, "ggHiNtuplizer", 5)'

# AllQCDPhoton30 , ggHiNtuplizerGED
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("AllQCDPhoton30", 1, "ggHiNtuplizerGED", 2)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("AllQCDPhoton30", 1, "ggHiNtuplizerGED", 3)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("AllQCDPhoton30", 1, "ggHiNtuplizerGED", 4)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("AllQCDPhoton30", 1, "ggHiNtuplizerGED", 5)'

# AllQCDPhoton30 , ggHiNtuplizer
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("EmEnrichedDijet30", 0, "ggHiNtuplizer", 2)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("EmEnrichedDijet30", 0, "ggHiNtuplizer", 3)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("EmEnrichedDijet30", 0, "ggHiNtuplizer", 4)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("EmEnrichedDijet30", 0, "ggHiNtuplizer", 5)'

# EmEnrichedDijet30 , ggHiNtuplizerGED
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("EmEnrichedDijet30", 0, "ggHiNtuplizerGED",2)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("EmEnrichedDijet30", 0, "ggHiNtuplizerGED",3)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("EmEnrichedDijet30", 0, "ggHiNtuplizerGED",4)'
root -l -b -q 'draw_isoDist_EPdep_sumIso.C+("EmEnrichedDijet30", 0, "ggHiNtuplizerGED",5)'

