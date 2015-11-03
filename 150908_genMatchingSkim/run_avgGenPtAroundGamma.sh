#void avgGenPtAroundGamma(const char* hiForestfileName="/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV.root",
#        TString outName="EmEnrichedDijet30", TString treePath="ggHiNtuplizer")


root -l -b -q 'avgGenPtAroundGamma.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV.root", "EmEnrichedDijet30", "ggHiNtuplizer")'
root -l -b -q 'avgGenPtAroundGamma.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV.root", "EmEnrichedDijet30", "ggHiNtuplizerGED")'
root -l -b -q 'avgGenPtAroundGamma.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV.root", "AllQCDPhoton30", "ggHiNtuplizer")'
root -l -b -q 'avgGenPtAroundGamma.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV.root", "AllQCDPhoton30", "ggHiNtuplizerGED")'
