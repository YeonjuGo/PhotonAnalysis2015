#void FilterEffi_data_1D(const char* fname="/afs/cern.ch/work/y/ygo/public/centrality/merged_centrality_MB_DATA750_RECO_150914.root", TString type="RECO_750", bool isAOD=0)

root -l -b -q 'FilterEffi_data_1D.C+("/u/user/goyeonju/files/centrality/Centrality_officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root", "HYDJET_5320" 0)'
root -l -b -q 'FilterEffi_data_1D.C+("/u/user/goyeonju/files/centrality/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611_CMSSW5320_byYJ.root", "2011_PbPb_5320",0)'
