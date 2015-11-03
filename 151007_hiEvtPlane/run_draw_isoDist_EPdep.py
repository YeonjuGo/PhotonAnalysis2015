import os

sample = ["AllQCDPhoton30","EmEnrichedDijet30"]
prompt = [0,1]
reco = ["ggHiNtuplizer","ggHiNtuplizerGED"]
var = ["pho_ecalClusterIsoR","pho_hcalRechitIsoR","pfcIso","pfnIso","pfpIso","pfcVsIso","pfnVsIso","pfpVsIso","sumIsoR","pfSumIso","pfSumVsIso"]
radius = [2,3,4,5]
evtNo = [15,21]

#evtNo = [2,8,15,21]


for c in reco:
    for e in evtNo:
        os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"AllQCDPhoton30\", 1, \"%s\", \"phoR9\",%d,0.0,1.0,50)'"%(c,e))
        os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"AllQCDPhoton30\", 1, \"%s\", \"phoHoverE\",%d,0.0,1.0,50)'"%(c,e))
        os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"AllQCDPhoton30\", 1, \"%s\", \"phoSigmaIEtaIEta\",%d,0.0,0.1,50)'"%(c,e))
        for d in radius:
            os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"AllQCDPhoton30\", 1, \"%s\", \"pho_trackIsoR%dPtCut20\",%d)'"%(c,d,e))
            for f in var:
                os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"AllQCDPhoton30\", 1, \"%s\", \"%s%d\",%d)'"%(c,f,d,e))


for c in reco:
    for e in evtNo:
        os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"EmEnrichedDijet30\", 0, \"%s\", \"phoR9\",%d,0.0,1.0,50)'"%(c,e))
        os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"EmEnrichedDijet30\", 0, \"%s\", \"phoHoverE\",%d,0.0,1.0,50)'"%(c,e))
        os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"EmEnrichedDijet30\", 0, \"%s\", \"phoSigmaIEtaIEta\",%d,0.0,0.1,50)'"%(c,e))
        for d in radius:
            os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"EmEnrichedDijet30\", 0, \"%s\", \"pho_trackIsoR%dPtCut20\",%d)'"%(c,d,e))
            for f in var:
                os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"EmEnrichedDijet30\", 0, \"%s\", \"%s%d\",%d)'"%(c,f,d,e))

        
"""
            os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"AllQCDPhoton30\", 1, \"%s\", \"pho_ecalClusterIsoR%d\",%d)'"%(j,i,k))
            os.system("root -l -b -q 'draw_isoDist_EPdep.C+(\"EmEnrichedDijet30\", 0, \"%s\", \"pho_ecalClusterIsoR%d\",%d)'"%(j,i,k))


"""
