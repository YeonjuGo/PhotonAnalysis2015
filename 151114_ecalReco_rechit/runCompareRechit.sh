 #void compareRechit(const char* f1="/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_global.root",
 #const char* f2="/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_multifit.root",
 #const char* cap="jet80Data")

#root -l -b -q 'compareRechit.C+("/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_global.root","/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_jet80_2011data_755p1_multifit.root", "jet80Data")'
#root -l -b -q 'compareRechit.C+("/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_AllQCDPhoton30_ecalGlobal_755p1_rechitPtMin0.root","/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_AllQCDPhoton30_ecalMultifit_755p1_rechitPtMin0.root", "photon30Data")'
#root -l -b -q 'compareRechit.C+("root://eoscms//eos/cms//store/group/phys_heavyions/mverweij/ECAL/FOREST/HiForestJet80Merged.root","/data/abaty/Pythia8_Dijet_pthat80_TuneCUETP8M1_5020GeV_HYDJET_MINBIAS_HiForest_PrivMC_MultiFit25ns/0.root", "jet80MC")'


#root -l 'compareRechit.C+("/afs/cern.ch/work/y/ygo/private/PRODUCTION/rechitChisqaure_test/CMSSW_7_5_5_patch1/src/EcalRecoStudy_chiSquare/test/HiForestAOD_jet80_multifit_chisquare.root","/afs/cern.ch/work/y/ygo/private/PRODUCTION/rechitChisqaure_test/CMSSW_7_5_5_patch1/src/EcalRecoStudy_chiSquare/test/HiForestAOD_photon_multifit_chisquare.root", "jet80_photon_data")'
root -l 'compareRechit.C+("/afs/cern.ch/work/y/ygo/private/PRODUCTION/rechitChisqaure_test/CMSSW_7_5_5_patch1/src/EcalRecoStudy_chiSquare/test/HiForestAOD_photon_multifit_chisquare.root","/afs/cern.ch/work/y/ygo/public/ecalLocalReco/photonMC_n300.root", "photon_data_mc")'
