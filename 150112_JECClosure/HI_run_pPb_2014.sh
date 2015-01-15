#in MIT
#The samples used to derive the forward pPb JEC have been re-forested using the derived JEC
: << 'bbb'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPblocalDBJECSampleClosureTest/HiForest_pPbJECClosure15.root","test/pPblocalDBJECSampleClosureTest15.root",30,84110.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPblocalDBJECSampleClosureTest/HiForest_pPbJECClosure30.root","test/pPblocalDBJECSampleClosureTest30.root",50,46533.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPblocalDBJECSampleClosureTest/HiForest_pPbJECClosure50.root","test/pPblocalDBJECSampleClosureTest50.root",80,13144.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPblocalDBJECSampleClosureTest/HiForest_pPbJECClosure80.root","test/pPblocalDBJECSampleClosureTest80.root",120,3096.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPblocalDBJECSampleClosureTest/HiForest_pPbJECClosure120.root","test/pPblocalDBJECSampleClosureTest120.root",9999,733.0)'

list=`echo ./test/pPblocalDB*.root`
hadd merged_pPblocalDBJECSampleClosureTest.root $list
bbb

#in KNU
#pPb MIX produced by Alex
#: << 'aaa'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root","test/MIX_pPb30.root",50,102400.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root","test/MIX_pPb50.root",80,39656.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root","test/MIX_pPb80.root",120,10157.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root","test/MIX_pPb120.root",170,2517.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton170.root","test/MIX_pPb170.root",9999,649.0)'

list=`echo ./test/MIX_pPb*.root`
hadd merged_MIX_pPb.root $list
#aaa

#in Korea University
#pPb 2013 old sample to test macro
:<<'ccc'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/PA2013_pyquen_allQCDPhoton30to50_forestv85.root","test/old_pPb30.root",50,56669.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/PA2013_pyquen_allQCDPhoton50to80_forestv85.root","test/old_pPb50.root",80,41906.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/PA2013_pyquen_allQCDPhoton80to120_forestv85.root","test/old_pPb80.root",120,12044.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/PA2013_pyquen_allQCDPhoton120to9999_forestv85.root","test/old_pPb120.root",9999,4481.0)'

#list=`echo ./test/old_pPb*.root`
#hadd merged_old_pPb.root $list
ccc

: << 'bbb'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/pA_Pyquen_allQCDPhoton30_hiForest2_53x_2013-18-14-1922.root","test/old_pPb30_53x_2013-18-14-1922.root",50,56669.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/pA_Pyquen_allQCDPhoton50_hiForest2_53x_2013-18-14-1922.root","test/old_pPb50_53x_2013-18-14-1922.root",80,41906.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/pA_Pyquen_allQCDPhoton80_hiForest2_53x_2013-18-14-1922.root","test/old_pPb80_53x_2013-18-14-1922.root",120,12044.0)'
root -l -b -q 'HI_weight_pPb_leadingPho.C+("/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/pA_Pyquen_allQCDPhoton120_hiForest2_53x_2013-18-14-1922.root","test/old_pPb120_53x_2013-18-14-1922.root",9999,4481.0)'

list=`echo ./test/old_pPb*_53x_2013-18-14-1922.root`
hadd merged_old_pPb_53x_2013-18-14-1922.root $list
bbb

