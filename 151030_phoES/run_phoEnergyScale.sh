#void phoEnergyScale_eta(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer", int     ptThr=20)
root -l -b -q 'phoEnergyScale_eta_centDep.C+("AllQCDPhoton30","ggHiNtuplizer")'
root -l -b -q 'phoEnergyScale_eta_centDep.C+("AllQCDPhoton30","ggHiNtuplizerGED")'

#void phoEnergyScale_dphi(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer", int ptThr=20, int evpOrder=2)
root -l -b -q 'phoEnergyScale_dphi_centDep.C+("AllQCDPhoton30","ggHiNtuplizer",20,2)'
root -l -b -q 'phoEnergyScale_dphi_centDep.C+("AllQCDPhoton30","ggHiNtuplizerGED",20,2)'

#void phoEnergyScale_pt(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer")
root -l -b -q 'phoEnergyScale_pt_centDep.C+("AllQCDPhoton30","ggHiNtuplizer")'
root -l -b -q 'phoEnergyScale_pt_centDep.C+("AllQCDPhoton30","ggHiNtuplizerGED")'

#void phoEnergyScale_r9_centDep(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer",int ptThr=20)
root -l -b -q 'phoEnergyScale_r9_centDep.C+("AllQCDPhoton30","ggHiNtuplizer")'
root -l -b -q 'phoEnergyScale_r9_centDep.C+("AllQCDPhoton30","ggHiNtuplizerGED")'

#void phoEnergyScale_brem_centDep(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer",int ptThr=20)
root -l -b -q 'phoEnergyScale_brem_centDep.C+("AllQCDPhoton30","ggHiNtuplizer")'
root -l -b -q 'phoEnergyScale_brem_centDep.C+("AllQCDPhoton30","ggHiNtuplizerGED")'


#################### SingleGammaFlatPt10To200
#void phoEnergyScale_eta(TString sample = "SingleGammaFlatPt10To200", TString treePath="ggHiNtuplizer", int ptThr=20)
root -l -b -q 'phoEnergyScale_eta.C+("SingleGammaFlatPt10To200","ggHiNtuplizer")'
root -l -b -q 'phoEnergyScale_eta.C+("SingleGammaFlatPt10To200","ggHiNtuplizerGED")'

#void phoEnergyScale_dphi(TString sample = "SingleGammaFlatPt10To200", TString treePath="ggHiNtuplizer", int ptThr=20, int evpOrder=2)
##### photon gun sample has no centrality info (just two photons per event)

#void phoEnergyScale_pt(TString sample = "SingleGammaFlatPt10To200", TString treePath="ggHiNtuplizer")
root -l -b -q 'phoEnergyScale_pt.C+("SingleGammaFlatPt10To200","ggHiNtuplizer")'
root -l -b -q 'phoEnergyScale_pt.C+("SingleGammaFlatPt10To200","ggHiNtuplizerGED")'

#void phoEnergyScale_r9(TString sample = "SingleGammaFlatPt10To200", TString treePath="ggHiNtuplizer",int ptThr=20)
root -l -b -q 'phoEnergyScale_r9.C+("SingleGammaFlatPt10To200","ggHiNtuplizer")'
root -l -b -q 'phoEnergyScale_r9.C+("SingleGammaFlatPt10To200","ggHiNtuplizerGED")'

#void phoEnergyScale_brem(TString sample = "SingleGammaFlatPt10To200", TString treePath="ggHiNtuplizer",int ptThr=20)
root -l -b -q 'phoEnergyScale_brem.C+("SingleGammaFlatPt10To200","ggHiNtuplizer")'
root -l -b -q 'phoEnergyScale_brem.C+("SingleGammaFlatPt10To200","ggHiNtuplizerGED")'


