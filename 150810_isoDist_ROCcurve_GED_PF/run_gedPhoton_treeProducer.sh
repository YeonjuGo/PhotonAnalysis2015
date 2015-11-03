#void gedPhotonAnalyzer(const char* hiForestfileName="HiForest.root", TString outName="AllQCDPhoton30", TString treePath="ggHiNtuplizer", float ptThr=30, condition cond_ = noC)
#enum condition {noC, hoeC, sigmaC, hoeAndSigmaC}

############ photon above 30
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 30)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 30)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 30)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 30)'

#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 30, hoeC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 30, hoeC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 30, hoeC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 30, hoeC, 100)'

#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 30, sigmaC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 30, sigmaC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 30, sigmaC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 30, sigmaC, 100)'
#
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 30, hoeAndSigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 30, hoeAndSigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 30, hoeAndSigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 30, hoeAndSigmaC, 80)'

############ photon above 40
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 40)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 40)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 40)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 40)'

#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 40, hoeC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 40, hoeC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 40, hoeC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 40, hoeC, 100)'
#
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 40, sigmaC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 40, sigmaC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 40, sigmaC, 100)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 40, sigmaC, 100)'
#
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 40, hoeAndSigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 40, hoeAndSigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 40, hoeAndSigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 40, hoeAndSigmaC, 80)'
#
############## photon above 60
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 60)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 60)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 60)'
root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 60)'

#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 60, hoeC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 60, hoeC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 60, hoeC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 60, hoeC, 80)'
#
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 60, sigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 60, sigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 60, sigmaC, 80)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 60, sigmaC, 80)'
#
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizer", 60, hoeAndSigmaC, 50)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_forest_pfisoAdded.root","AllQCDPhoton30", "ggHiNtuplizerGED", 60, hoeAndSigmaC, 50)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizer", 60, hoeAndSigmaC, 50)'
#root -l -b -q 'gedPhotonAnalyzer.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_forest_pfisoAdded.root","EmEnrichedDijet30", "ggHiNtuplizerGED", 60, hoeAndSigmaC, 50)'
