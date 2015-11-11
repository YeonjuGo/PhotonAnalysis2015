import os

sample = ["AllQCDPhoton30","EmEnrichedDijet30"]
prompt = [0,1]
reco = ["ggHiNtuplizer","ggHiNtuplizerGED"]
var = ["pho_ecalClusterIsoR","pho_hcalRechitIsoR","pfcIso","pfnIso","pfpIso","pfcVsIso","pfnVsIso","pfpVsIso","sumIso","pfSumIso","pfSumVsIso"]
radius = [2,3,4,5]
evtNo = [2,8,15,21]


for c in reco:
    for e in evtNo:
        for d in radius:
            os.system("root -l -b -q 'draw_isoDist_EPdep_sumIso_exceptPfIso.C+(\"AllQCDPhoton30\", 1, \"%s\",%d,%d)'"%(c,d,e))
            os.system("root -l -b -q 'draw_isoDist_EPdep_sumIso_exceptPfIso.C+(\"EmEnrichedDijet30\", 0, \"%s\",%d,%d)'"%(c,d,e))

