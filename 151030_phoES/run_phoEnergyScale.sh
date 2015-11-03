#void phoEnergyScale_eta(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer", int     ptThr=20)

root -l -b -q 'phoEnergyScale_eta.C+("AllQCDPhoton30","ggHiNtuplizer",20)'
root -l -b -q 'phoEnergyScale_eta.C+("AllQCDPhoton30","ggHiNtuplizer",40)'
root -l -b -q 'phoEnergyScale_eta.C+("AllQCDPhoton30","ggHiNtuplizer",60)'
#root -l -b -q 'phoEnergyScale_eta.C+("AllQCDPhoton30","ggHiNtuplizerGED",20)'
#root -l -b -q 'phoEnergyScale_eta.C+("AllQCDPhoton30","ggHiNtuplizerGED",40)'
#root -l -b -q 'phoEnergyScale_eta.C+("AllQCDPhoton30","ggHiNtuplizerGED",60)'

#void phoEnergyScale_dphi(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer", int ptThr=20, int evpOrder=2)

root -l -b -q 'phoEnergyScale_dphi.C+("AllQCDPhoton30","ggHiNtuplizer",20,2)'
#root -l -b -q 'phoEnergyScale_dphi.C+("AllQCDPhoton30","ggHiNtuplizerGED",20,2)'
#root -l -b -q 'phoEnergyScale_dphi.C+("AllQCDPhoton30","ggHiNtuplizer",20,3)'
#root -l -b -q 'phoEnergyScale_dphi.C+("AllQCDPhoton30","ggHiNtuplizerGED",20,3)'

#void phoEnergyScale_pt(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer")

root -l -b -q 'phoEnergyScale_pt.C+("AllQCDPhoton30","ggHiNtuplizer")'
#root -l -b -q 'phoEnergyScale_pt.C+("AllQCDPhoton30","ggHiNtuplizerGED")'

