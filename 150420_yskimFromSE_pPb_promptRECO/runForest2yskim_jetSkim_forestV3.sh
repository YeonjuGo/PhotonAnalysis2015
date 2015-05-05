#void forest2yskim_jetSkim_forestV3(TString inputFile_="forestFiles/pA/pA_photonSkimForest_v85_skimPhotonPt50_eta1.5.root",
#                                   std::string MinbiasFname = "skim_trackJet_minbiasTrackJet_mc.root",
#                                   float cutphotonPt  = 35,  // default value dropped to 35GeV  for later photon energy smearing/scaling
#                                   std::string outname = "testPhotonSkim.root",
#                                   sampleType colli=kPADATA,
#                                   bool doMix = false,
#                                   bool doJetResCorrection = 1,  // jet energy correction is done by default from Oct 19th (YS)
#				    int smearingCentBin = -1, //0=0-10%, 1=10-30%, 2=30-50%, 3=50-100%, 4=0-30%, 5=30-100%  : Jet pT and phi smearing!
#                                   float addFlatJetEnergyRes = 0,
#                                   bool useGenJetColl = 0
#                                   )



#enum sampleType { kHIDATA, kHIMC, kPPDATA, kPPMC, kPADATA, kPAMC};
#                       0    1      2       3       4        5



# # pbpb data
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pbpb/HiForestPhoton-v7-noDuplicate.root",       "yskimmedFiles/minbiasSkim_HIMinBias_merged_1_1000.root",    35,  "yskim_HiForestPhoton-v7-noDuplicate.root",  0, 1)'

# # pbpb MC
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pbpb/qcdAllPhoton30to50_genPhotonPtCut40_allCent.root",       "yskimmedFiles/minbiasSkim_HydjetDrum03_HiForest_v05_merged_test02.root",    35,  "yskim_qcdAllPhoton30to50_genPhotonPtCut40_allCent.root", 1, 1)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pbpb/qcdAllPhoton50to80_genPhotonPtCut40_allCent.root",       "yskimmedFiles/minbiasSkim_HydjetDrum03_HiForest_v05_merged_test02.root" ,   35,  "yskim_qcdAllPhoton50to80_genPhotonPtCut40_allCent.root", 1, 1)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pbpb/qcdAllPhoton80to9999_genPhotonPtCut40_allCent.root",       "yskimmedFiles/minbiasSkim_HydjetDrum03_HiForest_v05_merged_test02.root",    35,  "yskim_qcdAllPhoton80to9999_genPhotonPtCut40_allCent.root", 1, 1)'

# # pp data
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/pp_photonSKimForest_v85.root",       "no minbias", 35,   "yskimmedFiles/yskim_pp_photonSKimForest_v85_jetNoSmearing.root",2,0,1)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/pp_photonSKimForest_v85.root",       "no minbias", 35,   "yskimmedFiles/yskim_pp_photonSKimForest_v85_jetSmearingCent5099.root",2,0, 1, 3)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/pp_photonSKimForest_v85.root",       "no minbias", 35,   "yskimmedFiles/yskim_pp_photonSKimForest_v85_jetSmearingCent3050.root",2,0, 1, 2)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/pp_photonSKimForest_v85.root",       "no minbias", 35,   "yskimmedFiles/yskim_pp_photonSKimForest_v85_jetSmearingCent1030.root",2,0, 1, 1)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/pp_photonSKimForest_v85.root",       "no minbias", 35,   "yskimmedFiles/yskim_pp_photonSKimForest_v85_jetSmearingCent0010.root",2,0, 1, 0)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/pp_photonSKimForest_v85.root",       "no minbias", 35,   "yskimmedFiles/yskim_pp_photonSKimForest_v85_jetSmearingCent3099.root",2,0, 1, 5)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/pp_photonSKimForest_v85.root",       "no minbias", 35,   "yskimmedFiles/yskim_pp_photonSKimForest_v85_jetSmearingCent0030.root",2,0, 1, 4)'


# # pp MC. photon hiForest
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/merged_allQCDPhoton30to50_genPhotonPtCut30_CMSSW538HIp2.root",       "no min bias", 35,  "yskim_merged_allQCDPhoton30to50_genPhotonPtCut30_CMSSW538HIp2.root",3,0)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/merged_allQCDPhoton50to80_genPhotonPtCut30_CMSSW538HIp2.root",       "no min bias", 35,  "yskim_merged_allQCDPhoton50to80_genPhotonPtCut30_CMSSW538HIp2.root",3,0)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/merged_allQCDPhoton80to120_genPhotonPtCut30_CMSSW538HIp2.root",       "no min bias", 35,  "yskim_merged_allQCDPhoton80to120_genPhotonPtCut30_CMSSW538HIp2.root",3,0)'
# root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/merged_allQCDPhoton120to9999_genPhotonPtCut30_CMSSW538HIp2.root",       "no min bias", 35,  "yskim_merged_allQCDPhoton120to9999_genPhotonPtCut30_CMSSW538HIp2.root",3,0)'
# pp MC. background hiForest
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/emDijet30to50_CMSSW538HIp2.root",       "minbiasSkim_H02_merged_hiforest_v68_prod12_16HFBin.root", 35,  "yskim_emDijet30to50_CMSSW538HIp2.root",3,0)'
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/emDijet50to80_CMSSW538HIp2.root",       "minbiasSkim_H02_merged_hiforest_v68_prod12_16HFBin.root", 35,  "yskim_emDijet50to80_CMSSW538HIp2.root",3,0)'
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/emDijet80to120_CMSSW538HIp2.root",       "minbiasSkim_H02_merged_hiforest_v68_prod12_16HFBin.root", 35,  "yskim_emDijet80to120_CMSSW538HIp2.root",3,0)'
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("forestFiles/pp/emDijet120to9999_CMSSW538HIp2.root",       "minbiasSkim_H02_merged_hiforest_v68_prod12_16HFBin.root", 35,  "yskim_emDijet120to9999_CMSSW538HIp2.root",3,0)'


#  pA data, photon hiForest
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C++("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_DATA_photon30trig_localJEC_v1.root", "", 0, 1, 35,  "yskimmedFiles/yskim_pA_DATA_photonForest.root",4,0)'
# pA MC.  photon hiForest
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C++("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton30.root", "", 50,  102400, 35, "yskimmedFiles/yskim_pA_AllQCDPhoton3050.root", 5,0)'
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton50.root", "", 80, 39656, 35, "yskimmedFiles/yskim_pA_AllQCDPhoton5080.root", 5,0)'
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton80.root", "", 120, 10157, 35, "yskimmedFiles/yskim_pA_AllQCDPhoton80120.root", 5,0)'
root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton120.root", "", 170, 2517,35, "yskimmedFiles/yskim_pA_AllQCDPhoton120170.root", 5,0)'
#root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton170.root", "", 9999, 649, 35, "yskimmedFiles/yskim_pA_AllQCDPhoton1709999.root", 5,0)'
