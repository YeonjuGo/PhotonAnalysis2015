#void isoPhoEff_isoDep(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer",TString isoType = "sumIsoR", int rad=4, int evpOrder=2)

root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizer","sumIsoR",4,2)'
root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizerGED","sumIsoR",4,2)'
root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizer","pfSumVsIso",4,2)'
root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizerGED","pfSumVsIso",4,2)'

root -l -b -q 'isoPhoEff_isoDep.C+("EmEnrichedDijet30","ggHiNtuplizer","sumIsoR",4,2)'
root -l -b -q 'isoPhoEff_isoDep.C+("EmEnrichedDijet30","ggHiNtuplizerGED","sumIsoR",4,2)'
root -l -b -q 'isoPhoEff_isoDep.C+("EmEnrichedDijet30","ggHiNtuplizer","pfSumVsIso",4,2)'
root -l -b -q 'isoPhoEff_isoDep.C+("EmEnrichedDijet30","ggHiNtuplizerGED","pfSumVsIso",4,2)'



#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizer","sumIsoR",4,1)'
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizer","sumIsoR",4,3)'
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizer","sumIsoR",4,4)'
#
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizerGED","sumIsoR",4,1)'
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizerGED","sumIsoR",4,3)'
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizerGED","sumIsoR",4,4)'
#
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizer","pfSumVsIso",4,1)'
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizer","pfSumVsIso",4,3)'
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizer","pfSumVsIso",4,4)'
#
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizerGED","pfSumVsIso",4,1)'
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizerGED","pfSumVsIso",4,3)'
#root -l -b -q 'isoPhoEff_isoDep.C+("AllQCDPhoton30","ggHiNtuplizerGED","pfSumVsIso",4,4)'
