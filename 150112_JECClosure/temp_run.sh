#HI_weight_pPb.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat15_5020GeV.root" "test/pho_pPb15.root" 30 199048.0 
#HI_weight_pPb.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat30_5020GeV.root" "test/pho_pPb30.root" 50 12900.0
#HI_weight_pPb.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat50_5020GeV.root" "test/pho_pPb50.root" 80 1467.0
#HI_weight_pPb.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat80_5020GeV.root" "test/pho_pPb80.root" 120 145.0
#HI_weight_pPb.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat120_5020GeV.root" "test/pho_pPb120.root" 170 13.0
#HI_weight_pPb.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat170_5020GeV.root" "test/pho_pPb170.root" 9999 2.0

#HI_weight_pPb_leadingPho.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat15_5020GeV.root" "test/leadingpho_pPb15.root" 30 199048.0
#HI_weight_pPb_leadingPho.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat30_5020GeV.root" "test/leadingpho_pPb30.root" 50 12900.0
#HI_weight_pPb_leadingPho.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat50_5020GeV.root" "test/leadingpho_pPb50.root" 80 1467.0
#HI_weight_pPb_leadingPho.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat80_5020GeV.root" "test/leadingpho_pPb80.root" 120 145.0
#HI_weight_pPb_leadingPho.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat120_5020GeV.root" "test/leadingpho_pPb120.root" 170 13.0
#HI_weight_pPb_leadingPho.exe "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pp_pPbstyle_MC_5020GeV/pp_pPbstyle_MC_AllQCDPhotons_PtHat170_5020GeV.root" "test/leadingpho_pPb170.root" 9999 2.0

root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD30/JEC_pPbAllQCD30.root","test/pPb30_photonThreshold35_temp.root",9999,102400.0)'
#root -l -b -q 'HI_weight_pPb_leadingPho_changeThreshold.C++("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD30/JEC_pPbAllQCD30.root","test/leadingpho40_pPb30_photonThreshold35_temp.root",9999,102400.0)'
#root -l -b -q 'HI_weight_pPb_leadingPho.C++("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD30/JEC_pPbAllQCD30.root","test/leadingpho40_pPb30_photonThreshold25_temp.root",9999,46533.0)'
#root -l -b -q 'HI_weight_pPb_leadingPho.C++("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root","test/leadingpho40_pPb30_temp.root",9999,102400.0)'

: << 'sss'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root","test/leadingpho40_pPb50_bugFixed.root",80,39656.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root","test/leadingpho40_pPb80_bugFixed.root",120,10157.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root","test/leadingpho40_pPb120_bugFixed.root",170,2517.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton170.root","test/leadingpho40_pPb170_bugFixed.root",9999,649.0)'

list=`echo ./test/leadingpho40_pPb*_bugFixed.root`
hadd merged_pPb_leadingpho40_bugFixed.root $list
sss
