import os

out_f = open("test.sh",'w')

reco = ["ggHiNtuplizer","ggHiNtuplizerGED"]
var = ["pho_ecalClusterIsoR","pho_hcalRechitIsoR","pfcIso","pfnIso","pfpIso","pfcVsIso","pfnVsIso","pfpVsIso","sumIsoR","pfSumIso","pfSumVsIso"]
radius = [3]
#radius = [2,3,4,5]
evtNo = [2,8,15,21]

for c in reco:
    for e in evtNo:
        out_f.write("root -l -b -q 'draw_isoDist_EPdep_inclusivePho.C+(\"%s\", \"phoR9\",%d,0.0,1.0,50)'\n"%(c,e))
        out_f.write("root -l -b -q 'draw_isoDist_EPdep_inclusivePho.C+(\"%s\", \"phoHoverE\",%d,0.0,1.0,50)'\n"%(c,e))
        out_f.write("root -l -b -q 'draw_isoDist_EPdep_inclusivePho.C+(\"%s\", \"phoSigmaIEtaIEta\",%d,0.0,0.1,50)'\n"%(c,e))
        for d in radius:
            out_f.write("root -l -b -q 'draw_isoDist_EPdep_inclusivePho.C+(\"%s\", \"pho_trackIsoR%dPtCut20\",%d)'\n"%(c,d,e))
            for f in var:
                out_f.write("root -l -b -q 'draw_isoDist_EPdep_inclusivePho.C+(\"%s\", \"%s%d\",%d)'\n"%(c,f,d,e))
        
out_f.close()

