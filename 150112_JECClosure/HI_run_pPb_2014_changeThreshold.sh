
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD30/JEC_pPbAllQCD30.root","test/nonMIX_pPb30to9999_phoThr35.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD50/JEC_pPbAllQCD50.root","test/nonMIX_pPb50to9999_phoThr35.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD80/JEC_pPbAllQCD80.root","test/nonMIX_pPb80to9999_phoThr35.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD120/JEC_pPbAllQCD120.root","test/nonMIX_pPb120to9999_phoThr35.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD170/JEC_pPbAllQCD170.root","test/nonMIX_pPb170to9999_phoThr35.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_leadingPho_afterChangeThr.C++("test/nonMIX_pPb30to9999_phoThr35.root","test/leading_nonMIX_pPb30to9999_phoThr35.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_leadingPho_afterChangeThr.C++("test/nonMIX_pPb50to9999_phoThr35.root","test/leading_nonMIX_pPb50to9999_phoThr35.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_leadingPho_afterChangeThr.C++("test/nonMIX_pPb80to9999_phoThr35.root","test/leading_nonMIX_pPb80to9999_phoThr35.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_leadingPho_afterChangeThr.C++("test/nonMIX_pPb120to9999_phoThr35.root","test/leading_nonMIX_pPb120to9999_phoThr35.root",9999,1.0)'


: << 'sss'
# to check the HI_weight_pPb_changeThreshold.C macro
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root","test/MIX_pPb30to9999_phoThr35test.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root","test/MIX_pPb50to9999_phoThr35test.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root","test/MIX_pPb80to9999_phoThr35test.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root","test/MIX_pPb120to9999_phoThr35test.root",9999,1.0)'
root -l -b -q 'HI_weight_pPb_changeThreshold.C++("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton170.root","test/MIX_pPb170to9999_phoThr35test.root",9999,1.0)'

root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root","test/leadingpho40_pPb50_bugFixed.root",80,39656.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root","test/leadingpho40_pPb80_bugFixed.root",120,10157.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root","test/leadingpho40_pPb120_bugFixed.root",170,2517.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton170.root","test/leadingpho40_pPb170_bugFixed.root",9999,649.0)'

list=`echo ./test/leadingpho40_pPb*_bugFixed.root`
hadd merged_pPb_leadingpho40_bugFixed.root $list
sss
