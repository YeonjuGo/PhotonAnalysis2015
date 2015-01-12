#in KNU
#pPb MIX produced by Alex
: << 'END'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root","test/leadingpho40_pPb30_bugFixed.root",50,102400.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root","test/leadingpho40_pPb50_bugFixed.root",80,39656.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root","test/leadingpho40_pPb80_bugFixed.root",120,10157.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root","test/leadingpho40_pPb120_bugFixed.root",170,2517.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton170.root","test/leadingpho40_pPb170_bugFixed.root",9999,649.0)'

list=`echo ./test/leadingpho40_pPb*_bugFixed.root`
hadd merged_pPb_leadingpho40_bugFixed.root $list
END

#in Korea University
#pPb 2013 old sample to test macro

root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/pA_Pyquen_allQCDPhoton30_hiForest2_53x_2013-18-14-1922.root","test/old_pPb30.root",50,56669./50385)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/pA_Pyquen_allQCDPhoton50_hiForest2_53x_2013-18-14-1922.root","test/old_pPb50.root",80,41906./114136)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/pA_Pyquen_allQCDPhoton80_hiForest2_53x_2013-18-14-1922.root","test/old_pPb80.root",120,12044./103562)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/pA_Pyquen_allQCDPhoton120_hiForest2_53x_2013-18-14-1922.root","test/old_pPb120.root",170,4481./151511)'

list=`echo ./test/old_pPb*.root`
hadd merged_old_pPb.root $list


