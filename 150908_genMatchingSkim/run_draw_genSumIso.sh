#void avgGenPtAroundGamma(const char* hiForestfileName="/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV.root",
#        TString outName="EmEnrichedDijet30", TString treePath="ggHiNtuplizer")

root -l -b -q 'draw_genSumIso.C+("EmEnrichedDijet30", "ggHiNtuplizer",20,"genSumIso2>0&&abs(phoEta)<=1.44")'
root -l -b -q 'draw_genSumIso.C+("EmEnrichedDijet30", "ggHiNtuplizerGED",20,"genSumIso2>0&&abs(phoEta)<=1.44")'
root -l -b -q 'draw_genSumIso.C+("AllQCDPhoton30", "ggHiNtuplizer",20,"genSumIso2>0&&abs(phoEta)<=1.44")'
root -l -b -q 'draw_genSumIso.C+("AllQCDPhoton30", "ggHiNtuplizerGED",20,"genSumIso2>0&&abs(phoEta)<=1.44")'

root -l -b -q 'draw_genSumIso.C+("EmEnrichedDijet30", "ggHiNtuplizer",20,"genSumIso2>0&&abs(phoEta)>=1.566&&abs(phoEta)<2.0")'
root -l -b -q 'draw_genSumIso.C+("EmEnrichedDijet30", "ggHiNtuplizerGED",20,"genSumIso2>0&&abs(phoEta)>=1.566&&abs(phoEta)<2.0")'
root -l -b -q 'draw_genSumIso.C+("AllQCDPhoton30", "ggHiNtuplizer",20,"genSumIso2>0&&abs(phoEta)>=1.566&&abs(phoEta)<2.0")'
root -l -b -q 'draw_genSumIso.C+("AllQCDPhoton30", "ggHiNtuplizerGED",20,"genSumIso2>0&&abs(phoEta)>=1.566&&abs(phoEta)<2.0")'

root -l -b -q 'draw_genSumIso.C+("EmEnrichedDijet30", "ggHiNtuplizer",20,"genSumIso2>0&&abs(phoEta)>2.0")'
root -l -b -q 'draw_genSumIso.C+("EmEnrichedDijet30", "ggHiNtuplizerGED",20,"genSumIso2>0&&abs(phoEta)>2.0")'
root -l -b -q 'draw_genSumIso.C+("AllQCDPhoton30", "ggHiNtuplizer",20,"genSumIso2>0&&abs(phoEta)>2.0")'
root -l -b -q 'draw_genSumIso.C+("AllQCDPhoton30", "ggHiNtuplizerGED",20,"genSumIso2>0&&abs(phoEta)>2.0")'
