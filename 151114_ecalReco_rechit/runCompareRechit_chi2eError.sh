 #void compareRechit(const char* f1="/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_global.root",
 #const char* f2="/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_multifit.root",
 #const char* cap="jet80Data")

#root -l 'compareRechit_chi2eError.C+("","", "")'
#root -l 'compareRechit_chi2eError.C+("/afs/cern.ch/work/y/ygo/private/PRODUCTION/rechitChisqaure_test/CMSSW_7_5_5_patch1/src/EcalRecoStudy_chiSquare/test/HiForestAOD_jet80_multifit_chisquare.root","/afs/cern.ch/work/y/ygo/private/PRODUCTION/rechitChisqaure_test/CMSSW_7_5_5_patch1/src/EcalRecoStudy_chiSquare/test/HiForestAOD_photon_multifit_chisquare.root", "jet80_photon_data")'
#root -l 'compareRechit_chi2eError.C+("/afs/cern.ch/work/y/ygo/private/PRODUCTION/rechitChisqaure_test/CMSSW_7_5_5_patch1/src/EcalRecoStudy_chiSquare/test/HiForestAOD_photon_multifit_chisquare.root","/afs/cern.ch/work/y/ygo/public/ecalLocalReco/photonMC_n1000.root", "photon_data_mc")'

root -l 'compareRechit_chi2eError.C+("/afs/cern.ch/work/y/ygo/private/PRODUCTION/CMSSW_7_5_5_patch2/src/EcalReco/test/HiForestAOD_photonData_global_chisquare.root","/afs/cern.ch/work/y/ygo/private/PRODUCTION/CMSSW_7_5_5_patch2/src/EcalReco/test/HiForestAOD_photonData_multifit_chisquare.root", "photonData_global_multifit")'