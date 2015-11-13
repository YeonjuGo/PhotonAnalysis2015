#void FilterEffi_data_1D(const char* fname="/afs/cern.ch/work/y/ygo/public/centrality/merged_centrality_MB_DATA750_RECO_150914.root", TString type="RECO_750", bool isAOD=0)

root -l -b -q 'FilterEffi_data_1D.C+("/store/group/phys_heavyions/dgulhan/pp_minbiasSkim_forest_53x_2013-08-15-0155/pp_minbiasSkim_forest_53x_2013-08-15-0155.root", "ppDATA_MinBiasSkim_53x", 0, 0)'
