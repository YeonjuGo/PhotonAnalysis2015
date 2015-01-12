
root -l -b -q 'HI_weight_pPb_leadingPho_changeThreshold.C++("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root","test/leadingpho40_pPb30_photonThreshold35.root",50,102400.0)'
: << 'sss'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root","test/leadingpho40_pPb50_bugFixed.root",80,39656.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root","test/leadingpho40_pPb80_bugFixed.root",120,10157.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root","test/leadingpho40_pPb120_bugFixed.root",170,2517.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton170.root","test/leadingpho40_pPb170_bugFixed.root",9999,649.0)'

list=`echo ./test/leadingpho40_pPb*_bugFixed.root`
hadd merged_pPb_leadingpho40_bugFixed.root $list
sss
