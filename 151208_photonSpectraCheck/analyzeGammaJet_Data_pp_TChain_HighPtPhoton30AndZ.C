/*
 * macro to analyze and save di-Photon, di-Electron and di-Muon spectrum
 */

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TH1D.h>
#include <TList.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "../treeUtil.h"
#include "../histoUtil.h"
#include "../EventMatcher.h"

const int MAXJETS = 500;
const float cutPt = 10;
const float cutEta = 1.4791;
const float cutDeltaR = 0.15;
const long MAXTREESIZE = 200000000000; // set maximum tree size from 10 GB to 100 GB, so that the code does not switch to a new file after 10 GB7

const double eleMass = 0.000511;
const double muMass  = 0.105658;

const bool useHIPhotonIsolation = true;
const bool looseIsolation  = true;
#define PI 3.141592653589
const float awayRange = PI * 7./8.;

void analyzeGammaJet_Data_pp_TChain(const char* inputFile, const char* outputFile = "out_analyzeDiPhoEleMu_Data_HI.root");
Double_t getDR  ( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
Double_t getDETA(Double_t eta1, Double_t eta2);

void analyzeGammaJet_Data_pp_TChain(const char* inputFile, const char* outputFile)
{
       std::cout<<"running analyzeGammaJet_Data_pp_TChain()"<<std::endl;
       std::cout<<"inputFile   = "<< inputFile <<std::endl;
       std::cout<<"outputFile  = "<< outputFile <<std::endl;

       TChain* treeHLT   = new TChain("hltanalysis/HltTree");
       TChain* tree      = new TChain("ggHiNtuplizer/EventTree");
       TChain* treePho   = new TChain("ggHiNtuplizer/EventTree");
       TChain* treeEvent = new TChain("ggHiNtuplizer/EventTree");
       TChain* treeJet   = new TChain("ak4PFJetAnalyzer/t");

       int numFiles = 72;
       const char* fileNames[numFiles] =
       {
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_1.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_10.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_12.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_13.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_14.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_15.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_16.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_17.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_18.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_19.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_20.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_21.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_22.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_23.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_24.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_25.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_26.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_27.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_28.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_29.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_31.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_32.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_33.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_34.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_36.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_38.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_39.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_40.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_41.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_42.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_43.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_44.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_45.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_46.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_47.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_48.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_49.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_5.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_50.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_51.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_52.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_53.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_54.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_55.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_57.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_58.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_59.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_60.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_61.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_64.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_65.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_67.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_68.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_69.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_70.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_72.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_73.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_74.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_75.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_76.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_77.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_78.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_79.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_80.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_81.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_82.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_83.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_84.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_85.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_86.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_87.root",
               "/afs/cern.ch/user/k/katatar/eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtPhoton30AndZ/CRAB_UserFiles/crab_pp5Tev_HighPtPhoton30AndZ_PromptReco_take2/151201_211359/0000/HiForestAOD_9.root"
       };

       for (int i=0; i<numFiles; ++i)
       {
           treeHLT->Add(fileNames[i]);
           treePho->Add(fileNames[i]);
           treeEvent->Add(fileNames[i]);
           treeJet->Add(fileNames[i]);
           tree->Add(fileNames[i]);
       }

       treeHLT->SetBranchStatus("*",0);     // disable all branches
       treeHLT->SetBranchStatus("HLT_HI*SinglePhoton*Eta*v1*",1);     // enable photon branches
       treeHLT->SetBranchStatus("HLT_HI*DoublePhoton*Eta*v1*",1);     // enable photon branches
       treeHLT->SetBranchStatus("*DoubleMu*",1);     // enable muon branches

       treeEvent->SetBranchStatus("*",0);
       treeEvent->SetBranchStatus("run",1);
       treeEvent->SetBranchStatus("event",1);
       treeEvent->SetBranchStatus("lumis",1);

       treePho->SetBranchStatus("*",0);
       treePho->SetBranchStatus("run",1);
       treePho->SetBranchStatus("event",1);
       treePho->SetBranchStatus("lumis",1);
       treePho->SetBranchStatus("nPho",1);
       treePho->SetBranchStatus("pho*",1);

       treeJet->SetBranchStatus("*",0);        // disable all branches
       treeJet->SetBranchStatus("nref",1);     // enable jet branches
       treeJet->SetBranchStatus("jtpt",1);     // enable jet branches
       treeJet->SetBranchStatus("jteta",1);     // enable jet branches
       treeJet->SetBranchStatus("jtphi",1);     // enable jet branches

       Int_t nref;
       Float_t jtpt[MAXJETS];
       Float_t jteta[MAXJETS];
       Float_t jtphi[MAXJETS];

       treeJet->SetBranchAddress("nref",&nref);
       treeJet->SetBranchAddress("jtpt",jtpt);
       treeJet->SetBranchAddress("jteta",jteta);
       treeJet->SetBranchAddress("jtphi",jtphi);
       
       // event information
       Int_t run, lumis;
       Long64_t event;
       treeEvent->SetBranchAddress("run", &run);
       treeEvent->SetBranchAddress("event", &event);
       treeEvent->SetBranchAddress("lumis", &lumis);

       // RECO photons
       Int_t nPho;
       std::vector<float>* phoE=0;
       std::vector<float>* phoEt=0;
       std::vector<float>* phoEta=0;
       std::vector<float>* phoPhi=0;
       std::vector<float>* pho_ecalClusterIsoR4=0;
       std::vector<float>* pho_hcalRechitIsoR4=0;
       std::vector<float>* pho_trackIsoR4PtCut20=0;
       std::vector<float>* phoR9=0;
       std::vector<float>* phoHoverE=0;
       std::vector<float>* phoSigmaIEtaIEta=0;
       std::vector<float>* phoSigmaIEtaIEta_2012=0;
       std::vector<float>* phoSCRawE=0;
       std::vector<float>* phoE5x5=0;
       std::vector<float>* phoSCEtaWidth=0;
       std::vector<float>* phoSCPhiWidth=0;
       std::vector<float>* pho_swissCrx=0;
       std::vector<float>* pho_seedTime=0;
       std::vector<float>* pfcIso4=0;
       std::vector<float>* pfnIso4=0;

       // RECO electrons
       Int_t nEle;
       std::vector<int>*   eleCharge=0;
       std::vector<float>* eleEn=0;
       std::vector<float>* elePt=0;
       std::vector<float>* eleEta=0;
       std::vector<float>* elePhi=0;
       std::vector<float>* eleHoverE=0;
       std::vector<float>* eleSigmaIEtaIEta=0;
       std::vector<float>* eleSigmaIEtaIEta_2012=0;
       std::vector<float>* eleSigmaIPhiIPhi=0;
       std::vector<float>* eleEoverPInv=0;
       std::vector<float>* eledEtaAtVtx=0;
       std::vector<float>* eledPhiAtVtx=0;
       std::vector<float>* eleD0=0;
       std::vector<float>* eleDz=0;
       std::vector<float>* eleTrkPt=0;

       // RECO muons
       Int_t nMu;
       std::vector<int>*   muCharge=0;
       std::vector<float>* muPt=0;
       std::vector<float>* muEta=0;
       std::vector<float>* muPhi=0;
       std::vector<int>*   muIsGood=0;
       std::vector<float>* muD0=0;
       std::vector<float>* muDz=0;
       std::vector<float>* muChi2NDF=0;
       std::vector<float>* muInnerD0=0;
       std::vector<float>* muInnerDz=0;
       std::vector<float>* muIsoTrk=0;
       std::vector<float>* muPFChIso=0;
       std::vector<float>* muPFPhoIso=0;
       std::vector<float>* muPFNeuIso=0;
       std::vector<float>* muPFPUIso=0;

       tree->SetBranchStatus("*",0);     // disable all branches
       tree->SetBranchStatus("nPho",1);     // enable photon branches
       tree->SetBranchStatus("pho*",1);     // enable photon branches
       tree->SetBranchAddress("nPho",&nPho);
       tree->SetBranchAddress("phoE",&phoE);
       tree->SetBranchAddress("phoEt",&phoEt);
       tree->SetBranchAddress("phoEta",&phoEta);
       tree->SetBranchAddress("phoPhi",&phoPhi);
       tree->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
       tree->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
       tree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
       tree->SetBranchAddress("phoR9",&phoR9);
       tree->SetBranchAddress("phoHoverE",&phoHoverE);
       tree->SetBranchAddress("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
       tree->SetBranchAddress("phoSigmaIEtaIEta_2012",&phoSigmaIEtaIEta_2012);
       tree->SetBranchAddress("phoSCRawE",&phoSCRawE);
       tree->SetBranchAddress("phoE5x5",&phoE5x5);
       tree->SetBranchAddress("phoSCEtaWidth",&phoSCEtaWidth);
       tree->SetBranchAddress("phoSCPhiWidth",&phoSCPhiWidth);
       tree->SetBranchAddress("pho_swissCrx",&pho_swissCrx);
       tree->SetBranchAddress("pho_seedTime",&pho_seedTime);
       tree->SetBranchAddress("pfcIso4",&pfcIso4);
       tree->SetBranchAddress("pfnIso4",&pfnIso4);

       tree->SetBranchStatus("nEle",1);     // enable electron branches
       tree->SetBranchStatus("ele*",1);     // enable electron branches
       tree->SetBranchAddress("nEle",&nEle);
       tree->SetBranchAddress("eleEn",&eleEn);
       tree->SetBranchAddress("eleCharge",&eleCharge);
       tree->SetBranchAddress("elePt",&elePt);
       tree->SetBranchAddress("eleEta",&eleEta);
       tree->SetBranchAddress("elePhi",&elePhi);
       tree->SetBranchAddress("eleHoverE",&eleHoverE);
       tree->SetBranchAddress("eleSigmaIEtaIEta",&eleSigmaIEtaIEta);
       tree->SetBranchAddress("eleSigmaIEtaIEta_2012",&eleSigmaIEtaIEta_2012);
       tree->SetBranchAddress("eleSigmaIPhiIPhi",&eleSigmaIPhiIPhi);
       tree->SetBranchAddress("eleEoverPInv",&eleEoverPInv);
       tree->SetBranchAddress("eledEtaAtVtx",&eledEtaAtVtx);
       tree->SetBranchAddress("eledPhiAtVtx",&eledPhiAtVtx);
       tree->SetBranchAddress("eleD0",&eleD0);
       tree->SetBranchAddress("eleDz",&eleDz);
       tree->SetBranchAddress("eleTrkPt",&eleTrkPt);

       tree->SetBranchStatus("nMu",1);     // enable muon branches
       tree->SetBranchStatus("mu*",1);     // enable muon branches
       tree->SetBranchAddress("nMu",&nMu);
       tree->SetBranchAddress("muCharge",&muCharge);
       tree->SetBranchAddress("muPt",&muPt);
       tree->SetBranchAddress("muEta",&muEta);
       tree->SetBranchAddress("muPhi",&muPhi);
       tree->SetBranchAddress("muIsGood",&muIsGood);
       tree->SetBranchAddress("muD0",&muD0);
       tree->SetBranchAddress("muDz",&muDz);
       tree->SetBranchAddress("muChi2NDF",&muChi2NDF);
       tree->SetBranchAddress("muInnerD0",&muInnerD0);
       tree->SetBranchAddress("muInnerDz",&muInnerDz);
       tree->SetBranchAddress("muIsoTrk",&muIsoTrk);
       tree->SetBranchAddress("muPFChIso",&muPFChIso);
       tree->SetBranchAddress("muPFPhoIso",&muPFPhoIso);
       tree->SetBranchAddress("muPFNeuIso",&muPFNeuIso);
       tree->SetBranchAddress("muPFPUIso",&muPFPUIso);

       TFile* output = new TFile(outputFile,"RECREATE");
       // output tree variables

       TTree *outputTreeHLT    = treeHLT->CloneTree(0);
       outputTreeHLT->SetName("hltTree");
       TTree *outputTreePho    = treePho->CloneTree(0);
       outputTreePho->SetName("photons");
       outputTreePho->SetTitle("Event data + photons");
       TTree *outputTreeJet    = treeJet->CloneTree(0);
       outputTreeJet->SetName("jets");
       TTree *outputTreeGammaJet = new TTree("gammaJet","leading photon-jet correlations");
       TTree *outputTreediPho  = new TTree("diphoton","");
       TTree *outputTreediEle  = new TTree("dielectron","");
       TTree *outputTreediMu   = new TTree("dimuon","");

       outputTreeHLT->SetMaxTreeSize(MAXTREESIZE);
       outputTreediPho->SetMaxTreeSize(MAXTREESIZE);
       outputTreediEle->SetMaxTreeSize(MAXTREESIZE);
       outputTreediMu->SetMaxTreeSize(MAXTREESIZE);

       Int_t phoIdx;
       std::vector<int>   jetIdx;
       std::vector<float> xjg;
       std::vector<float> deta;
       std::vector<float> dphi;
       std::vector<float> dR;
       std::vector<int>   insideJet;
       Int_t nJetin7over8;

       outputTreeGammaJet->Branch("phoIdx",&phoIdx);
       outputTreeGammaJet->Branch("jetIdx",&jetIdx);
       outputTreeGammaJet->Branch("xjg",&xjg);
       outputTreeGammaJet->Branch("deta",&deta);
       outputTreeGammaJet->Branch("dphi",&dphi);
       outputTreeGammaJet->Branch("dR",&dR);
       outputTreeGammaJet->Branch("insideJet",&insideJet);
       outputTreeGammaJet->Branch("nJetin7over8",&nJetin7over8);

       // output tree variables
       std::vector<float> phoEt_1;
       std::vector<float> phoEta_1;
       std::vector<float> phoPhi_1;
       std::vector<float> pho_ecalClusterIsoR4_1;
       std::vector<float> pho_hcalRechitIsoR4_1;
       std::vector<float> pho_trackIsoR4PtCut20_1;
       std::vector<float> phoR9_1;
       std::vector<float> phoHoverE_1;
       std::vector<float> phoSigmaIEtaIEta_1;
       std::vector<float> phoSigmaIEtaIEta_2012_1;
       std::vector<float> phoSCRawE_1;
       std::vector<float> phoE5x5_1;
       std::vector<float> phoSCEtaWidth_1;
       std::vector<float> phoSCPhiWidth_1;
       std::vector<float> pho_swissCrx_1;
       std::vector<float> pho_seedTime_1;
       std::vector<int>   matched_eleCharge_1;
       std::vector<int>   matched_eleIndex_1;
       std::vector<float> matched_elePt_1;
       std::vector<float> matched_eleEta_1;
       std::vector<float> matched_elePhi_1;
       std::vector<float> matched_eleTrkPt_1;

       std::vector<float> phoEt_2;
       std::vector<float> phoEta_2;
       std::vector<float> phoPhi_2;
       std::vector<float> pho_ecalClusterIsoR4_2;
       std::vector<float> pho_hcalRechitIsoR4_2;
       std::vector<float> pho_trackIsoR4PtCut20_2;
       std::vector<float> phoR9_2;
       std::vector<float> phoHoverE_2;
       std::vector<float> phoSigmaIEtaIEta_2;
       std::vector<float> phoSigmaIEtaIEta_2012_2;
       std::vector<float> phoSCRawE_2;
       std::vector<float> phoE5x5_2;
       std::vector<float> phoSCEtaWidth_2;
       std::vector<float> phoSCPhiWidth_2;
       std::vector<float> pho_swissCrx_2;
       std::vector<float> pho_seedTime_2;
       std::vector<int>   matched_eleCharge_2;
       std::vector<int>   matched_eleIndex_2;
       std::vector<float> matched_elePt_2;
       std::vector<float> matched_eleEta_2;
       std::vector<float> matched_elePhi_2;
       std::vector<float> matched_eleTrkPt_2;

       std::vector<float> vSum_gg_M;
       std::vector<float> vSum_gg_Energy;
       std::vector<float> vSum_gg_Pt;
       std::vector<float> vSum_gg_Eta;
       std::vector<float> vSum_gg_Phi;

       std::vector<float> vSum_gg_M_elePt;
       std::vector<float> vSum_gg_M_eleTrkPt;

       outputTreediPho->Branch("nPho",&nPho);
       outputTreediPho->Branch("phoEt_1",&phoEt_1);
       outputTreediPho->Branch("phoEta_1",&phoEta_1);
       outputTreediPho->Branch("phoPhi_1",&phoPhi_1);
       outputTreediPho->Branch("pho_ecalClusterIsoR4_1",&pho_ecalClusterIsoR4_1);
       outputTreediPho->Branch("pho_hcalRechitIsoR4_1",&pho_hcalRechitIsoR4_1);
       outputTreediPho->Branch("pho_trackIsoR4PtCut20_1",&pho_trackIsoR4PtCut20_1);
       outputTreediPho->Branch("phoR9_1",&phoR9_1);
       outputTreediPho->Branch("phoHoverE_1",&phoHoverE_1);
       outputTreediPho->Branch("phoSigmaIEtaIEta_1",&phoSigmaIEtaIEta_1);
       outputTreediPho->Branch("phoSigmaIEtaIEta_2012_1",&phoSigmaIEtaIEta_2012_1);
       outputTreediPho->Branch("phoSCRawE_1",&phoSCRawE_1);
       outputTreediPho->Branch("phoE5x5_1",&phoE5x5_1);
       outputTreediPho->Branch("phoSCEtaWidth_1",&phoSCEtaWidth_1);
       outputTreediPho->Branch("phoSCPhiWidth_1",&phoSCPhiWidth_1);
       outputTreediPho->Branch("pho_swissCrx_1",&pho_swissCrx_1);
       outputTreediPho->Branch("pho_seedTime_1",&pho_seedTime_1);
       outputTreediPho->Branch("matched_eleCharge_1",&matched_eleCharge_1);
       outputTreediPho->Branch("matched_eleIndex_1",&matched_eleIndex_1);
       outputTreediPho->Branch("matched_elePt_1",&matched_elePt_1);
       outputTreediPho->Branch("matched_eleEta_1",&matched_eleEta_1);
       outputTreediPho->Branch("matched_elePhi_1",&matched_elePhi_1);
       outputTreediPho->Branch("matched_eleTrkPt_1",&matched_eleTrkPt_1);

       outputTreediPho->Branch("phoEt_2",&phoEt_2);
       outputTreediPho->Branch("phoEta_2",&phoEta_2);
       outputTreediPho->Branch("phoPhi_2",&phoPhi_2);
       outputTreediPho->Branch("pho_ecalClusterIsoR4_2",&pho_ecalClusterIsoR4_2);
       outputTreediPho->Branch("pho_hcalRechitIsoR4_2",&pho_hcalRechitIsoR4_2);
       outputTreediPho->Branch("pho_trackIsoR4PtCut20_2",&pho_trackIsoR4PtCut20_2);
       outputTreediPho->Branch("phoR9_2",&phoR9_2);
       outputTreediPho->Branch("phoHoverE_2",&phoHoverE_2);
       outputTreediPho->Branch("phoSigmaIEtaIEta_2",&phoSigmaIEtaIEta_2);
       outputTreediPho->Branch("phoSigmaIEtaIEta_2012_2",&phoSigmaIEtaIEta_2012_2);
       outputTreediPho->Branch("phoSCRawE_2",&phoSCRawE_2);
       outputTreediPho->Branch("phoE5x5_2",&phoE5x5_2);
       outputTreediPho->Branch("phoSCEtaWidth_2",&phoSCEtaWidth_2);
       outputTreediPho->Branch("phoSCPhiWidth_2",&phoSCPhiWidth_2);
       outputTreediPho->Branch("pho_swissCrx_2",&pho_swissCrx_2);
       outputTreediPho->Branch("pho_seedTime_2",&pho_seedTime_2);
       outputTreediPho->Branch("matched_eleCharge_2",&matched_eleCharge_2);
       outputTreediPho->Branch("matched_eleIndex_2",&matched_eleIndex_2);
       outputTreediPho->Branch("matched_elePt_2",&matched_elePt_2);
       outputTreediPho->Branch("matched_eleEta_2",&matched_eleEta_2);
       outputTreediPho->Branch("matched_elePhi_2",&matched_elePhi_2);
       outputTreediPho->Branch("matched_eleTrkPt_2",&matched_eleTrkPt_2);

       outputTreediPho->Branch("vSum_M",&vSum_gg_M);
       outputTreediPho->Branch("vSum_Energy",&vSum_gg_Energy);
       outputTreediPho->Branch("vSum_Pt",&vSum_gg_Pt);
       outputTreediPho->Branch("vSum_Eta",&vSum_gg_Eta);
       outputTreediPho->Branch("vSum_Phi",&vSum_gg_Phi);

       outputTreediPho->Branch("vSum_M_elePt",&vSum_gg_M_elePt);
       outputTreediPho->Branch("vSum_M_eleTrkPt",&vSum_gg_M_eleTrkPt);

       std::vector<int>   eleCharge_1;
       std::vector<float> elePt_1;
       std::vector<float> eleEta_1;
       std::vector<float> elePhi_1;
       std::vector<float> eleHoverE_1;
       std::vector<float> eleSigmaIEtaIEta_1;
       std::vector<float> eleSigmaIEtaIEta_2012_1;
       std::vector<float> eleSigmaIPhiIPhi_1;
       std::vector<float> eleEoverPInv_1;
       std::vector<float> eledEtaAtVtx_1;
       std::vector<float> eledPhiAtVtx_1;
       std::vector<float> eleD0_1;
       std::vector<float> eleDz_1;

       std::vector<int>   eleCharge_2;
       std::vector<float> elePt_2;
       std::vector<float> eleEta_2;
       std::vector<float> elePhi_2;
       std::vector<float> eleHoverE_2;
       std::vector<float> eleSigmaIEtaIEta_2;
       std::vector<float> eleSigmaIEtaIEta_2012_2;
       std::vector<float> eleSigmaIPhiIPhi_2;
       std::vector<float> eleEoverPInv_2;
       std::vector<float> eledEtaAtVtx_2;
       std::vector<float> eledPhiAtVtx_2;
       std::vector<float> eleD0_2;
       std::vector<float> eleDz_2;

       std::vector<float> vSum_ee_M;
       std::vector<float> vSum_ee_Energy;
       std::vector<float> vSum_ee_Pt;
       std::vector<float> vSum_ee_Eta;
       std::vector<float> vSum_ee_Phi;

       outputTreediEle->Branch("nEle",&nEle);
       outputTreediEle->Branch("eleCharge_1",&eleCharge_1);
       outputTreediEle->Branch("elePt_1",&elePt_1);
       outputTreediEle->Branch("eleEta_1",&eleEta_1);
       outputTreediEle->Branch("elePhi_1",&elePhi_1);
       outputTreediEle->Branch("eleHoverE_1",&eleHoverE_1);
       outputTreediEle->Branch("eleSigmaIEtaIEta_1",&eleSigmaIEtaIEta_1);
       outputTreediEle->Branch("eleSigmaIEtaIEta_2012_1",&eleSigmaIEtaIEta_2012_1);
       outputTreediEle->Branch("eleSigmaIPhiIPhi_1",&eleSigmaIPhiIPhi_1);
       outputTreediEle->Branch("eleEoverPInv_1",&eleEoverPInv_1);
       outputTreediEle->Branch("eledEtaAtVtx_1",&eledEtaAtVtx_1);
       outputTreediEle->Branch("eledPhiAtVtx_1",&eledPhiAtVtx_1);
       outputTreediEle->Branch("eleD0_1",&eleD0_1);
       outputTreediEle->Branch("eleDz_1",&eleDz_1);

       outputTreediEle->Branch("eleCharge_2",&eleCharge_2);
       outputTreediEle->Branch("elePt_2",&elePt_2);
       outputTreediEle->Branch("eleEta_2",&eleEta_2);
       outputTreediEle->Branch("elePhi_2",&elePhi_2);
       outputTreediEle->Branch("eleHoverE_2",&eleHoverE_2);
       outputTreediEle->Branch("eleSigmaIEtaIEta_2",&eleSigmaIEtaIEta_2);
       outputTreediEle->Branch("eleSigmaIEtaIEta_2012_2",&eleSigmaIEtaIEta_2012_2);
       outputTreediEle->Branch("eleSigmaIPhiIPhi_2",&eleSigmaIPhiIPhi_2);
       outputTreediEle->Branch("eleEoverPInv_2",&eleEoverPInv_2);
       outputTreediEle->Branch("eledEtaAtVtx_2",&eledEtaAtVtx_2);
       outputTreediEle->Branch("eledPhiAtVtx_2",&eledPhiAtVtx_2);
       outputTreediEle->Branch("eleD0_2",&eleD0_2);
       outputTreediEle->Branch("eleDz_2",&eleDz_2);

       outputTreediEle->Branch("vSum_M",&vSum_ee_M);
       outputTreediEle->Branch("vSum_Energy",&vSum_ee_Energy);
       outputTreediEle->Branch("vSum_Pt",&vSum_ee_Pt);
       outputTreediEle->Branch("vSum_Eta",&vSum_ee_Eta);
       outputTreediEle->Branch("vSum_Phi",&vSum_ee_Phi);

       // RECO muons
       std::vector<int>   muCharge_1;
       std::vector<float> muPt_1;
       std::vector<float> muEta_1;
       std::vector<float> muPhi_1;
       std::vector<int>   muIsGood_1;
       std::vector<float> muD0_1;
       std::vector<float> muDz_1;
       std::vector<float> muChi2NDF_1;
       std::vector<float> muInnerD0_1;
       std::vector<float> muInnerDz_1;
       std::vector<float> muIsoTrk_1;
       std::vector<float> muPFChIso_1;
       std::vector<float> muPFPhoIso_1;
       std::vector<float> muPFNeuIso_1;
       std::vector<float> muPFPUIso_1;

       std::vector<int>   muCharge_2;
       std::vector<float> muPt_2;
       std::vector<float> muEta_2;
       std::vector<float> muPhi_2;
       std::vector<int>   muIsGood_2;
       std::vector<float> muD0_2;
       std::vector<float> muDz_2;
       std::vector<float> muChi2NDF_2;
       std::vector<float> muInnerD0_2;
       std::vector<float> muInnerDz_2;
       std::vector<float> muIsoTrk_2;
       std::vector<float> muPFChIso_2;
       std::vector<float> muPFPhoIso_2;
       std::vector<float> muPFNeuIso_2;
       std::vector<float> muPFPUIso_2;

       std::vector<float> vSum_mumu_M;
       std::vector<float> vSum_mumu_Energy;
       std::vector<float> vSum_mumu_Pt;
       std::vector<float> vSum_mumu_Eta;
       std::vector<float> vSum_mumu_Phi;

       outputTreediMu->Branch("nMu",&nMu);
       outputTreediMu->Branch("muCharge_1",&muCharge_1);
       outputTreediMu->Branch("muPt_1",&muPt_1);
       outputTreediMu->Branch("muEta_1",&muEta_1);
       outputTreediMu->Branch("muPhi_1",&muPhi_1);
       outputTreediMu->Branch("muIsGood_1",&muIsGood_1);
       outputTreediMu->Branch("muD0_1",&muD0_1);
       outputTreediMu->Branch("muDz_1",&muDz_1);
       outputTreediMu->Branch("muChi2NDF_1",&muChi2NDF_1);
       outputTreediMu->Branch("muInnerD0_1",&muInnerD0_1);
       outputTreediMu->Branch("muInnerDz_1",&muInnerDz_1);
       outputTreediMu->Branch("muIsoTrk_1",&muIsoTrk_1);
       outputTreediMu->Branch("muPFChIso_1",&muPFChIso_1);
       outputTreediMu->Branch("muPFPhoIso_1",&muPFPhoIso_1);
       outputTreediMu->Branch("muPFNeuIso_1",&muPFNeuIso_1);
       outputTreediMu->Branch("muPFPUIso_1",&muPFPUIso_1);

       outputTreediMu->Branch("muCharge_2",&muCharge_2);
       outputTreediMu->Branch("muPt_2",&muPt_2);
       outputTreediMu->Branch("muEta_2",&muEta_2);
       outputTreediMu->Branch("muPhi_2",&muPhi_2);
       outputTreediMu->Branch("muIsGood_2",&muIsGood_2);
       outputTreediMu->Branch("muD0_2",&muD0_2);
       outputTreediMu->Branch("muDz_2",&muDz_2);
       outputTreediMu->Branch("muChi2NDF_2",&muChi2NDF_2);
       outputTreediMu->Branch("muInnerD0_2",&muInnerD0_2);
       outputTreediMu->Branch("muInnerDz_2",&muInnerDz_2);
       outputTreediMu->Branch("muIsoTrk_2",&muIsoTrk_2);
       outputTreediMu->Branch("muPFChIso_2",&muPFChIso_2);
       outputTreediMu->Branch("muPFPhoIso_2",&muPFPhoIso_2);
       outputTreediMu->Branch("muPFNeuIso_2",&muPFNeuIso_2);
       outputTreediMu->Branch("muPFPUIso_2",&muPFPUIso_2);

       outputTreediMu->Branch("vSum_M",&vSum_mumu_M);
       outputTreediMu->Branch("vSum_Energy",&vSum_mumu_Energy);
       outputTreediMu->Branch("vSum_Pt",&vSum_mumu_Pt);
       outputTreediMu->Branch("vSum_Eta",&vSum_mumu_Eta);
       outputTreediMu->Branch("vSum_Phi",&vSum_mumu_Phi);

       EventMatcher* em = new EventMatcher();
       Long64_t duplicateEntries = 0;

       Long64_t entries = treeEvent->GetEntries();
       Long64_t entriesAnalyzed = 0;
       std::cout << "entries         = " << entries << std::endl;
       std::cout<< "Loop : ggHiNtuplizer/EventTree" <<std::endl;
       for (Long64_t j_entry=0; j_entry<entries; ++j_entry)
       {
           if (j_entry % 2000 == 0)  {
             std::cout << "current entry = " <<j_entry<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j_entry/entries*100<<" %"<<std::endl;
           }

           treeHLT->GetEntry(j_entry);
           treeEvent->GetEntry(j_entry);
           treePho->GetEntry(j_entry);
           treeJet->GetEntry(j_entry);
           tree->GetEntry(j_entry);

           bool eventAdded = em->addEvent(run,lumis,event,j_entry);
           if(!eventAdded) // this event is duplicate, skip this one.
           {
               duplicateEntries++;
               continue;
           }

           // photon-jet block
           phoIdx = -1;
           double maxPhoEt = -1;
           for(int i=0; i<nPho; ++i)
           {
               bool passedEtaCut = TMath::Abs(phoEta->at(i)) < cutEta ;

               bool passedSpikeRejection;
               bool passedIsolation;
               if (useHIPhotonIsolation) {
                   passedSpikeRejection = (phoSigmaIEtaIEta->at(i)  > 0.002  &&
                                           pho_swissCrx->at(i)      < 0.9    &&
                                           TMath::Abs(pho_seedTime->at(i) ) < 3);
                   if (looseIsolation){
                       passedIsolation = (pho_ecalClusterIsoR4->at(i)  < 4.2 &&
                                          pho_hcalRechitIsoR4->at(i)   < 2.2 &&
                                          pho_trackIsoR4PtCut20->at(i) < 2   &&
                                          phoHoverE->at(i) < 0.1);
                   }
                   else
                   {
                       passedIsolation = ((pho_ecalClusterIsoR4->at(i) + pho_hcalRechitIsoR4->at(i) + pho_trackIsoR4PtCut20->at(i)) < 1 &&
                               phoHoverE->at(i) < 0.1) ;
                   }
               }
               else // particle flow isolation
               {
                   // HI Photon Isolation may not be in pp HiForest.
                   // use particle flow isolation in that case.

                   passedSpikeRejection = (phoSigmaIEtaIEta->at(i) > 0.002 );
                   if (looseIsolation){
                       passedIsolation = (pfcIso4->at(i)   < 5 &&
                                          pfnIso4->at(i)   < 5 &&
                                          phoHoverE->at(i) < 0.1);
                   }
                   else
                   {
                       passedIsolation = ((pfcIso4->at(i) + pfnIso4->at(i)) < 20 && phoHoverE->at(i) < 0.1) ;
                   }
               }

               bool passedPurity = (phoSigmaIEtaIEta->at(i) < 0.01);

               if (!(passedEtaCut && passedSpikeRejection && passedIsolation && passedPurity)) continue;

               if (phoEt->at(i) > maxPhoEt)
               {
                   maxPhoEt = phoEt->at(i);
                   phoIdx = i;
               }
           }
           if (phoIdx == -1) continue;
           // photon-jet correlation)
           nJetin7over8 = 0;

           jetIdx.clear();
           xjg.clear();
           deta.clear();
           dphi.clear();
           dR.clear();
           insideJet.clear();
           for (int i=0; i<nref; ++i)
           {
               // cuts on jets will be applied during plotting

               float tmp_deta = getDETA(phoEta->at(phoIdx), jteta[i]);
               float tmp_dphi = getDPHI(phoPhi->at(phoIdx), jtphi[i]);
               if (tmp_dphi > awayRange)
                   nJetin7over8++;
               float tmp_dR   = getDR(phoEta->at(phoIdx), phoPhi->at(phoIdx), jteta[i], jtphi[i]);

               int tmp_insideJet;
               if(tmp_dR < 0.4) tmp_insideJet = 1;
               else             tmp_insideJet = 0;

               jetIdx.push_back(i);
               xjg.push_back((float)jtpt[i]/phoEt->at(phoIdx));
               deta.push_back(tmp_deta);
               dphi.push_back(tmp_dphi);
               dR.push_back(tmp_dR);
               insideJet.push_back(tmp_insideJet);
           }

           outputTreeHLT->Fill();
           outputTreePho->Fill();
           outputTreeJet->Fill();
           
           outputTreeGammaJet->Fill();
           outputTreediPho->Fill();
           outputTreediEle->Fill();
           outputTreediMu->Fill();
       }
       std::cout<<  "Loop ENDED : ggHiNtuplizer/EventTree" <<std::endl;
       std::cout << "entries            = " << entries << std::endl;
       std::cout << "duplicateEntries   = " << duplicateEntries << std::endl;
       std::cout << "entriesAnalyzed    = " << entriesAnalyzed << std::endl;
       std::cout << "outputTreeHLT->GetEntries()   = " << outputTreeHLT->GetEntries() << std::endl;
       std::cout << "outputTreePho->GetEntries()      = " << outputTreePho->GetEntries() << std::endl;
       std::cout << "outputTreeJet->GetEntries()   = " << outputTreeJet->GetEntries() << std::endl;
    
       std::cout << "outputTreeGammaJet->GetEntries() = " << outputTreeGammaJet->GetEntries() << std::endl;
       std::cout << "outputTreediPho->GetEntries()    = " << outputTreediPho->GetEntries() << std::endl;
       std::cout << "outputTreediEle->GetEntries()    = " << outputTreediEle->GetEntries() << std::endl;
       std::cout << "outputTreediMu->GetEntries()     = " << outputTreediMu->GetEntries() << std::endl;
  
       output->Write();
       output->Close();
}

int main(int argc, char** argv)
{
    if (argc == 3)    analyzeGammaJet_Data_pp_TChain(argv[1], argv[2]);
    if (argc == 2)    analyzeGammaJet_Data_pp_TChain(argv[1]);
    if (argc == 1)    return 1;
    return 0;
}


Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
  Double_t theDphi = getDPHI( phi1, phi2);
  Double_t theDeta = eta1 - eta2;
  return TMath::Sqrt ( theDphi*theDphi + theDeta*theDeta);
}

Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;

  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 )
    dphi = dphi + 2. * 3.141592653589;

  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
  }

  return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2){
    return eta1 - eta2;
}
