#void gedPhoton_effi_figures(TString treePath=="ggHiNtuplizer", float ptThr=30, condition cond_ = noC )

# pt above 30 GeV
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 30, noC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 30, noC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 30, hoeC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 30, hoeC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 30, sigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 30, sigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 30, hoeAndSigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 30, hoeAndSigmaC)'

# pt above 40 GeV
root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 40, noC)'
root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 40, noC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 40, hoeC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 40, hoeC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 40, sigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 40, sigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 40, hoeAndSigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 40, hoeAndSigmaC)'

# pt above 60 GeV
root -l -b -q 'gedPhoton_effi_figures.C++("ggHiNtuplizer", 60, noC)'
root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 60, noC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 60, hoeC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 60, hoeC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 60, sigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 60, sigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizer", 60, hoeAndSigmaC)'
#root -l -b -q 'gedPhoton_effi_figures.C+("ggHiNtuplizerGED", 60, hoeAndSigmaC)'


