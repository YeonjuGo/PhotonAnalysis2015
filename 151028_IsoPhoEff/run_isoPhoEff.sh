#void isoPhoEff(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer",TString isoType = "sumIsoR", float isoCut=1.0, int rad=4, int evpOrder=2)

root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizer","sumIsoR",1.0,4,1)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizer","sumIsoR",1.0,4,2)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizer","sumIsoR",1.0,4,3)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizer","sumIsoR",1.0,4,4)'

root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizerGED","sumIsoR",1.0,4,1)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizerGED","sumIsoR",1.0,4,2)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizerGED","sumIsoR",1.0,4,3)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizerGED","sumIsoR",1.0,4,4)'

root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizer","pfSumVsIso",3.0,4,1)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizer","pfSumVsIso",3.0,4,2)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizer","pfSumVsIso",3.0,4,3)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizer","pfSumVsIso",3.0,4,4)'

root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizerGED","pfSumVsIso",3.0,4,1)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizerGED","pfSumVsIso",3.0,4,2)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizerGED","pfSumVsIso",3.0,4,3)'
root -l -b -q 'isoPhoEff.C+("AllQCDPhoton30","ggHiNtuplizerGED","pfSumVsIso",3.0,4,4)'
