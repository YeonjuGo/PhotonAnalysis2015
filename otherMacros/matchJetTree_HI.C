#include "EventMatchingCMS.h"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <iostream>
#include <TStyle.h>
#include "TH1F.h"
#include "TMath.h"

//const TString AnaFilename = "HiForest_Ncoll_Dijet_pthat80_740pre8_MCHI2_74_V3_merged_forest_0.root";
const TString AnaFilename = "HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v5.root";
//const TString AnaFilename = "PyquenUnquenched_DiJet_pt30_PbPb_5020GeV_actuallyEmbedded_HiForest.root";
//const TString AnaFilename = "AllQCDPhoton30_PhotonFilter20GeV_eta3_TuneZ2_PbPb_5020GeV_actuallyEmbedded_HiForest.root";
const TString Ana3CaloTreename = "akPu3CaloJetAnalyzer/t";
const TString Ana4CaloTreename = "akPu4CaloJetAnalyzer/t";
const TString AnaPhotonTreename = "multiPhotonAnalyzer/photon";
const TString AnaHITreename = "hiEvtAnalyzer/HiTree";
const TString AnaHLTTreename = "hltanalysis/HltTree";
const TString AnaSkimTreename = "skimanalysis/HltTree";
const TString AnaTrkTreename = "anaTrack/trackTree";
const TString AnaGenTreename = "HiGenParticleAna/hi";

//const TString HLTFilename = "openHLT_20150508_HIJet80502_740F.root";
const TString HLTFilename = "openHLT_20150508_HIMinBias502_740F.root";

const int nBins = 200;
const double maxpt = 200;

void matchJetTree()
{

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  TFile *HLTFile = TFile::Open(HLTFilename);
  //  TTree *HLTTree = (TTree*)HLTFile->Get("HltTree");
  TTree *HLTTree = (TTree*)HLTFile->Get("hltbitanalysis/HltTree");

  ULong64_t hlt_event;
  Int_t hlt_run, hlt_lumi;

  const Int_t n2Trig = 5;
  const Int_t n3Trig = 5;
  const Int_t n4Trig = 5;
  const Int_t n5Trig = 5;
  const Int_t nPhotonTrig = 5;

  Float_t L1IsolEmEt[4];
  Float_t L1NIsolEmEt[4];
  Float_t L1CenJetEt[4];
  Float_t L1ForJetEt[4];
  Int_t L1_SingleMu7;
  Int_t L1_SingleEG12;
  Int_t L1_SingleEG20;

  Int_t HLT_PuAK2CaloJet40_v1;
  Int_t HLT_PuAK2CaloJet60_v1;
  Int_t HLT_PuAK2CaloJet80_v1;
  Int_t HLT_PuAK2CaloJet100_v1;
  Int_t HLT_PuAK2CaloJet120_v1;

  Int_t HLT_PuAK2CaloJet40_MID_v1;
  Int_t HLT_PuAK2CaloJet60_MID_v1;
  Int_t HLT_PuAK2CaloJet80_MID_v1;
  Int_t HLT_PuAK2CaloJet100_MID_v1;
  Int_t HLT_PuAK2CaloJet120_MID_v1;

  Int_t HLT_PuAK3CaloJet40_v1;
  Int_t HLT_PuAK3CaloJet60_v1;
  Int_t HLT_PuAK3CaloJet80_v1;
  Int_t HLT_PuAK3CaloJet100_v1;
  Int_t HLT_PuAK3CaloJet120_v1;

  Int_t HLT_PuAK4CaloJet40_v1;
  Int_t HLT_PuAK4CaloJet60_v1;
  Int_t HLT_PuAK4CaloJet80_v1;
  Int_t HLT_PuAK4CaloJet100_v1;
  Int_t HLT_PuAK4CaloJet120_v1;

  Int_t HLT_PuAK4CaloJet40_MID_v1;
  Int_t HLT_PuAK4CaloJet60_MID_v1;
  Int_t HLT_PuAK4CaloJet80_MID_v1;
  Int_t HLT_PuAK4CaloJet100_MID_v1;
  Int_t HLT_PuAK4CaloJet120_MID_v1;

  Int_t HLT_PuAK4CaloJet80_45_45_MID_v1;

  Int_t HLT_PuAK4CaloJet80_MID2_v1;
  Int_t HLT_PuAK4CaloJet80_2MID2_v1;
  Int_t HLT_PuAK4CaloJet80_MID3_v1;
  Int_t HLT_PuAK4CaloJet80_2MID3_v1;
  Int_t HLT_PuAK4CaloJet100_MID2_v1;
  Int_t HLT_PuAK4CaloJet100_2MID2_v1;
  Int_t HLT_PuAK4CaloJet100_MID3_v1;
  Int_t HLT_PuAK4CaloJet100_2MID3_v1;

  Int_t HLT_PuAK5CaloJet40_v1;
  Int_t HLT_PuAK5CaloJet60_v1;
  Int_t HLT_PuAK5CaloJet80_v1;
  Int_t HLT_PuAK5CaloJet100_v1;
  Int_t HLT_PuAK5CaloJet120_v1;

  Int_t HLT_HISinglePhoton10_v1;
  Int_t HLT_HISinglePhoton15_v1;
  Int_t HLT_HISinglePhoton20_v1;
  Int_t HLT_HISinglePhoton40_v1;
  Int_t HLT_HISinglePhoton60_v1;

  Int_t HLT_HISinglePhoton10_MID_v1;
  Int_t HLT_HISinglePhoton15_MID_v1;
  Int_t HLT_HISinglePhoton20_MID_v1;
  Int_t HLT_HISinglePhoton40_MID_v1;
  Int_t HLT_HISinglePhoton60_MID_v1;

  Int_t HLT_HIDoublePhoton10_v1;
  Int_t HLT_HIDoublePhoton15_v1;
  Int_t HLT_HIDoublePhoton20_v1;
  Int_t HLT_HIDoublePhoton30_v1;

  Int_t HLT_HIPhoton15_Photon10_v1;
  Int_t HLT_HIPhoton20_Photon10_v1;
  Int_t HLT_HIPhoton20_Photon15_v1;

  Int_t nTrueJet40 = 0;
  Int_t nTrueJet80 = 0;
  Int_t nTrueJet90 = 0;
  Int_t nTrueJet100 = 0;
  Int_t nTrueJet110 = 0;
  Int_t nTrueJet120 = 0;

  Int_t nRecoJet40 = 0;
  Int_t nRecoJet80 = 0;
  //  Int_t nTruePhoton10 = 0;

  Int_t n2CaloJet40 = 0;
  Int_t n2CaloJet60 = 0;
  Int_t n2CaloJet80 = 0;
  Int_t n2CaloJet100 = 0;
  Int_t n2CaloJet120 = 0;

  Int_t n2CaloJet40_MID = 0;
  Int_t n2CaloJet60_MID = 0;
  Int_t n2CaloJet80_MID = 0;
  Int_t n2CaloJet100_MID = 0;
  Int_t n2CaloJet120_MID = 0;

  Int_t n3CaloJet40 = 0;
  Int_t n3CaloJet60 = 0;
  Int_t n3CaloJet80 = 0;
  Int_t n3CaloJet100 = 0;
  Int_t n3CaloJet120 = 0;

  Int_t n4CaloJet40 = 0;
  Int_t n4CaloJet60 = 0;
  Int_t n4CaloJet80 = 0;
  Int_t n4CaloJet100 = 0;
  Int_t n4CaloJet120 = 0;

  Int_t n4CaloJet40_MID = 0;
  Int_t n4CaloJet60_MID = 0;
  Int_t n4CaloJet80_MID = 0;
  Int_t n4CaloJet100_MID = 0;
  Int_t n4CaloJet120_MID = 0;

  Int_t n4CaloJet80_45_45_MID = 0;
  Int_t n4CaloJet80_45_45_MID_no100 = 0;
  Int_t n4CaloJet80_45_45_MID_no120 = 0;

  Int_t n4CaloJet80_MID2 = 0;
  Int_t n4CaloJet80_2MID2 = 0;
  Int_t n4CaloJet80_MID3 = 0;
  Int_t n4CaloJet80_2MID3 = 0;
  Int_t n4CaloJet100_MID2 = 0;
  Int_t n4CaloJet100_2MID2 = 0;
  Int_t n4CaloJet100_MID3 = 0;
  Int_t n4CaloJet100_2MID3 = 0;

  Int_t n4CaloJet80_MID2_no100 = 0;
  Int_t n4CaloJet80_2MID2_no100 = 0;
  Int_t n4CaloJet80_MID3_no100 = 0;
  Int_t n4CaloJet80_2MID3_no100 = 0;
  Int_t n4CaloJet100_MID2_no120 = 0;
  Int_t n4CaloJet100_2MID2_no120 = 0;
  Int_t n4CaloJet100_MID3_no120 = 0;
  Int_t n4CaloJet100_2MID3_no120 = 0;

  Int_t n4CaloJet100_SinglePhoton40 = 0;
  Int_t n4CaloJet100_SinglePhoton40_MID = 0;

  Int_t n5CaloJet40 = 0;
  Int_t n5CaloJet60 = 0;
  Int_t n5CaloJet80 = 0;
  Int_t n5CaloJet100 = 0;
  Int_t n5CaloJet120 = 0;

  Int_t nPhoton10 = 0;
  Int_t nPhoton15 = 0;
  Int_t nPhoton20 = 0;
  Int_t nPhoton40 = 0;
  Int_t nPhoton60 = 0;

  Int_t nPhoton10_MID = 0;
  Int_t nPhoton15_MID = 0;
  Int_t nPhoton20_MID = 0;
  Int_t nPhoton40_MID = 0;
  Int_t nPhoton60_MID = 0;

  Int_t nDoublePhoton10 = 0;
  Int_t nDoublePhoton15 = 0;
  Int_t nDoublePhoton20 = 0;
  Int_t nDoublePhoton30 = 0;

  Int_t nPhoton15_10 = 0;
  Int_t nPhoton20_10 = 0;
  Int_t nPhoton20_15 = 0;

  Int_t nFake2CaloJet40 = 0;
  Int_t nFake2CaloJet60 = 0;
  Int_t nFake2CaloJet80 = 0;
  Int_t nFake2CaloJet100 = 0;
  Int_t nFake2CaloJet120 = 0;

  Int_t nFake2CaloJet40_MID = 0;
  Int_t nFake2CaloJet60_MID = 0;
  Int_t nFake2CaloJet80_MID = 0;
  Int_t nFake2CaloJet100_MID = 0;
  Int_t nFake2CaloJet120_MID = 0;

  Int_t nFake3CaloJet40 = 0;
  Int_t nFake3CaloJet60 = 0;
  Int_t nFake3CaloJet80 = 0;
  Int_t nFake3CaloJet100 = 0;
  Int_t nFake3CaloJet120 = 0;

  Int_t nFake4CaloJet40 = 0;
  Int_t nFake4CaloJet60 = 0;
  Int_t nFake4CaloJet80 = 0;
  Int_t nFake4CaloJet100 = 0;
  Int_t nFake4CaloJet120 = 0;

  Float_t maxMissedPt80 = -1;
  Float_t maxMissedPt100 = -1;
  Float_t maxMissedPt120 = -1;

  Float_t maxMissedEta80 = -1;
  Float_t maxMissedEta100 = -1;
  Float_t maxMissedEta120 = -1;

  Float_t maxMissedGenPt80 = -1;
  Float_t maxMissedGenPt100 = -1;
  Float_t maxMissedGenPt120 = -1;

  Float_t maxMissedGenEta80 = -1;
  Float_t maxMissedGenEta100 = -1;
  Float_t maxMissedGenEta120 = -1;

  TH1F* totPt80_h = new TH1F("totPt80_h", "totPt80_h", 30, 80, 110);
  TH1F* totPt100_h = new TH1F("totPt100_h", "totPt100_h", 20, 100, 120);
  TH1F* totPt120_h = new TH1F("totPt120_h", "totPt120_h", 20, 120, 140);

  TH1F* totIntPt80_h = new TH1F("totIntPt80_h", "totIntPt80_h", 30, 80, 110);
  TH1F* totIntPt100_h = new TH1F("totIntPt100_h", "totIntPt100_h", 20, 100, 120);
  TH1F* totIntPt120_h = new TH1F("totIntPt120_h", "totIntPt120_h", 20, 120, 140);

  TH1F* totEta80_h = new TH1F("totEta80_h", "totEta80_h", 40, -2.0, 2.0);
  TH1F* totEta100_h = new TH1F("totEta100_h", "totEta100_h", 40, -2.0, 2.0);
  TH1F* totEta120_h = new TH1F("totEta120_h", "totEta120_h", 40, -2.0, 2.0);

  TH1F* totHiBin80_h = new TH1F("totHiBin80_h", "totHiBin80_h", 100, -0.5, 199.5);
  TH1F* totHiBin100_h = new TH1F("totHiBin100_h", "totHiBin100_h", 100, -0.5, 199.5);
  TH1F* totHiBin120_h = new TH1F("totHiBin120_h", "totHiBin120_h", 100, -0.5, 199.5);

  TH1F* maxMissedPt80_h = new TH1F("maxMissedPt80_h", "maxMissedPt80_h", 30, 80, 110);
  TH1F* maxMissedPt100_h = new TH1F("maxMissedPt100_h", "maxMissedPt100_h", 20, 100, 120);
  TH1F* maxMissedPt120_h = new TH1F("maxMissedPt120_h", "maxMissedPt120_h", 20, 120, 140);

  TH1F* maxMissedIntPt80_h = new TH1F("maxMissedIntPt80_h", "maxMissedIntPt80_h", 30, 80, 110);
  TH1F* maxMissedIntPt100_h = new TH1F("maxMissedIntPt100_h", "maxMissedIntPt100_h", 20, 100, 120);
  TH1F* maxMissedIntPt120_h = new TH1F("maxMissedIntPt120_h", "maxMissedIntPt120_h", 20, 120, 140);

  TH1F* maxMissedEta80_h = new TH1F("maxMissedEta80_h", "maxMissedEta80_h", 40, -2.0, 2.0);
  TH1F* maxMissedEta100_h = new TH1F("maxMissedEta100_h", "maxMissedEta100_h", 40, -2.0, 2.0);
  TH1F* maxMissedEta120_h = new TH1F("maxMissedEta120_h", "maxMissedEta120_h", 40, -2.0, 2.0);

  TH1F* maxMissedHiBin80_h = new TH1F("maxMissedHiBin80_h", "maxMissedHiBin80_h", 100, -0.5, 199.5);
  TH1F* maxMissedHiBin100_h = new TH1F("maxMissedHiBin100_h", "maxMissedHiBin100_h", 100, -0.5, 199.5);
  TH1F* maxMissedHiBin120_h = new TH1F("maxMissedHiBin120_h", "maxMissedHiBin120_h", 100, -0.5, 199.5);


  Int_t nFake4CaloJet80_45_45_MID = 0;
  Int_t nFake4CaloJet80_45_45_MID_no100 = 0;
  Int_t nFake4CaloJet80_45_45_MID_no120 = 0;

  Int_t nFake4CaloJet40_MID = 0;
  Int_t nFake4CaloJet60_MID = 0;
  Int_t nFake4CaloJet80_MID = 0;
  Int_t nFake4CaloJet100_MID = 0;
  Int_t nFake4CaloJet120_MID = 0;

  Int_t nFake5CaloJet40 = 0;
  Int_t nFake5CaloJet60 = 0;
  Int_t nFake5CaloJet80 = 0;
  Int_t nFake5CaloJet100 = 0;
  Int_t nFake5CaloJet120 = 0;

  Int_t nFakePhoton10 = 0;
  Int_t nFakePhoton15 = 0;
  Int_t nFakePhoton20 = 0;
  Int_t nFakePhoton40 = 0;
  Int_t nFakePhoton60 = 0;

  Int_t nFakePhoton10_MID = 0;
  Int_t nFakePhoton15_MID = 0;
  Int_t nFakePhoton20_MID = 0;
  Int_t nFakePhoton40_MID = 0;
  Int_t nFakePhoton60_MID = 0;

  Int_t nFakeDoublePhoton10 = 0;
  Int_t nFakeDoublePhoton15 = 0;
  Int_t nFakeDoublePhoton20 = 0;
  Int_t nFakeDoublePhoton30 = 0;

  Int_t nFakePhoton15_10 = 0;
  Int_t nFakePhoton20_10 = 0;
  Int_t nFakePhoton20_15 = 0;

  TString ak2CaloName[n2Trig] = {"HLT_PuAK2CaloJet40_v1", "HLT_PuAK2CaloJet60_v1", "HLT_PuAK2CaloJet80_v1", "HLT_PuAK2CaloJet100_v1", "HLT_PuAK2CaloJet120_v1"};
  TString ak3CaloName[n3Trig] = {"HLT_PuAK3CaloJet40_v1", "HLT_PuAK3CaloJet60_v1", "HLT_PuAK3CaloJet80_v1", "HLT_PuAK3CaloJet100_v1", "HLT_PuAK3CaloJet120_v1"};
  TString ak4CaloName[n4Trig] = {"HLT_PuAK4CaloJet40_v1", "HLT_PuAK4CaloJet60_v1", "HLT_PuAK4CaloJet80_v1", "HLT_PuAK4CaloJet100_v1", "HLT_PuAK4CaloJet120_v1"};
  TString ak5CaloName[n5Trig] = {"HLT_PuAK5CaloJet40_v1", "HLT_PuAK5CaloJet60_v1", "HLT_PuAK5CaloJet80_v1", "HLT_PuAK5CaloJet100_v1", "HLT_PuAK5CaloJet120_v1"};

  TString ak2CaloName_MID[n2Trig] = {"HLT_PuAK2CaloJet40_MID_v1", "HLT_PuAK2CaloJet60_MID_v1", "HLT_PuAK2CaloJet80_MID_v1", "HLT_PuAK2CaloJet100_MID_v1", "HLT_PuAK2CaloJet120_MID_v1"};
  TString ak4CaloName_MID[n4Trig] = {"HLT_PuAK4CaloJet40_MID_v1", "HLT_PuAK4CaloJet60_MID_v1", "HLT_PuAK4CaloJet80_MID_v1", "HLT_PuAK4CaloJet100_MID_v1", "HLT_PuAK4CaloJet120_MID_v1"};

  TString ak2CaloTrkName[n2Trig] = {"HLT_PuAK2CaloTrkJet40_v1", "HLT_PuAK2CaloTrkJet60_v1", "HLT_PuAK2CaloTrkJet80_v1", "HLT_PuAK2CaloTrkJet100_v1", "HLT_PuAK2CaloTrkJet120_v1"};
  TString ak3CaloTrkName[n3Trig] = {"HLT_PuAK3CaloTrkJet40_v1", "HLT_PuAK3CaloTrkJet60_v1", "HLT_PuAK3CaloTrkJet80_v1", "HLT_PuAK3CaloTrkJet100_v1", "HLT_PuAK3CaloTrkJet120_v1"};
  TString ak4CaloTrkName[n4Trig] = {"HLT_PuAK4CaloTrkJet40_v1", "HLT_PuAK4CaloTrkJet60_v1", "HLT_PuAK4CaloTrkJet80_v1", "HLT_PuAK4CaloTrkJet100_v1", "HLT_PuAK4CaloTrkJet120_v1"};
  TString ak5CaloTrkName[n5Trig] = {"HLT_PuAK5CaloTrkJet40_v1", "HLT_PuAK5CaloTrkJet60_v1", "HLT_PuAK5CaloTrkJet80_v1", "HLT_PuAK5CaloTrkJet100_v1", "HLT_PuAK5CaloTrkJet120_v1"};

  TString ak2CaloTrkName_MID[n2Trig] = {"HLT_PuAK2CaloTrkJet40_MID_v1", "HLT_PuAK2CaloTrkJet60_MID_v1", "HLT_PuAK2CaloTrkJet80_MID_v1", "HLT_PuAK2CaloTrkJet100_MID_v1", "HLT_PuAK2CaloTrkJet120_MID_v1"};
  TString ak4CaloTrkName_MID[n4Trig] = {"HLT_PuAK4CaloTrkJet40_MID_v1", "HLT_PuAK4CaloTrkJet60_MID_v1", "HLT_PuAK4CaloTrkJet80_MID_v1", "HLT_PuAK4CaloTrkJet100_MID_v1", "HLT_PuAK4CaloTrkJet120_MID_v1"};

  TString ak4CaloName_ULTRAMID[8] = {"HLT_PuAK4CaloJet80_MID2_v1", "HLT_PuAK4CaloJet80_2MID2_v1", "HLT_PuAK4CaloJet80_MID3_v1", "HLT_PuAK4CaloJet80_2MID3_v1", "HLT_PuAK4CaloJet100_MID2_v1", "HLT_PuAK4CaloJet100_2MID2_v1", "HLT_PuAK4CaloJet100_MID3_v1", "HLT_PuAK4CaloJet100_2MID3_v1"};

  TString ak4CaloName_3Jet[1] = {"HLT_PuAK4CaloJet80_45_45_MID_v1"};

  TString ak2CaloGenName[n2Trig] = {"HLT_PuAK2CaloGenJet40_v1", "HLT_PuAK2CaloGenJet60_v1", "HLT_PuAK2CaloGenJet80_v1", "HLT_PuAK2CaloGenJet100_v1", "HLT_PuAK2CaloGenJet120_v1"};
  TString ak3CaloGenName[n3Trig] = {"HLT_PuAK3CaloGenJet40_v1", "HLT_PuAK3CaloGenJet60_v1", "HLT_PuAK3CaloGenJet80_v1", "HLT_PuAK3CaloGenJet100_v1", "HLT_PuAK3CaloGenJet120_v1"};
  TString ak4CaloGenName[n4Trig] = {"HLT_PuAK4CaloGenJet40_v1", "HLT_PuAK4CaloGenJet60_v1", "HLT_PuAK4CaloGenJet80_v1", "HLT_PuAK4CaloGenJet100_v1", "HLT_PuAK4CaloGenJet120_v1"};
  TString ak5CaloGenName[n5Trig] = {"HLT_PuAK5CaloGenJet40_v1", "HLT_PuAK5CaloGenJet60_v1", "HLT_PuAK5CaloGenJet80_v1", "HLT_PuAK5CaloGenJet100_v1", "HLT_PuAK5CaloGenJet120_v1"};

  TString photonName[nPhotonTrig] = {"HLT_HISinglePhoton10_v1", "HLT_HISinglePhoton15_v1", "HLT_HISinglePhoton20_v1", "HLT_HISinglePhoton40_v1", "HLT_HISinglePhoton60_v1"};
  TString photonName_MID[nPhotonTrig] = {"HLT_HISinglePhoton10_MID_v1", "HLT_HISinglePhoton15_MID_v1", "HLT_HISinglePhoton20_MID_v1", "HLT_HISinglePhoton40_MID_v1", "HLT_HISinglePhoton60_MID_v1"};

  TString photonDoubleName[nPhotonTrig+2] = {"HLT_HIDoublePhoton10_v1", "HLT_HIDoublePhoton15_v1", "HLT_HIDoublePhoton20_v1", "HLT_HIDoublePhoton30_v1", "HLT_HIPhoton15_Photon10_v1", "HLT_HIPhoton20_Photon10_v1", "HLT_HIPhoton20_Photon15_v1"}; 

  const std::string genLeadJt3Name[3] = {"genLeadJet3Pt_h", "genLeadJet3Eta_h", "genLeadJet3Phi_h"};
  const std::string recoLeadJt3Name[3] = {"recoLeadJet3Pt_h", "recoLeadJet3Eta_h", "recoLeadJet3Phi_h"};
  const std::string genLeadJt4Name[3] = {"genLeadJet4Pt_h", "genLeadJet4Eta_h", "genLeadJet4Phi_h"};
  const std::string recoLeadJt4Name[3] = {"recoLeadJet4Pt_h", "recoLeadJet4Eta_h", "recoLeadJet4Phi_h"};

  const Int_t bins[3] = {200, 50, 50};
  const Float_t lower[3] = {0.0, -5.0, (Float_t)-TMath::Pi()};
  const Float_t upper[3] = {200.0, 5.0, (Float_t)TMath::Pi()};

  TH1F* genLeadJet3_p[3];
  TH1F* recoLeadJet3_p[3];
  TH1F* genLeadJet4_p[3];
  TH1F* recoLeadJet4_p[3];

  for(Int_t iter = 0; iter < 3; iter++){
    genLeadJet3_p[iter] = new TH1F(genLeadJt3Name[iter].c_str(), genLeadJt3Name[iter].c_str(), bins[iter], lower[iter], upper[iter]);
    recoLeadJet3_p[iter] = new TH1F(recoLeadJt3Name[iter].c_str(), recoLeadJt3Name[iter].c_str(), bins[iter], lower[iter], upper[iter]);

    genLeadJet4_p[iter] = new TH1F(genLeadJt4Name[iter].c_str(), genLeadJt4Name[iter].c_str(), bins[iter], lower[iter], upper[iter]);
    recoLeadJet4_p[iter] = new TH1F(recoLeadJt4Name[iter].c_str(), recoLeadJt4Name[iter].c_str(), bins[iter], lower[iter], upper[iter]);
  }


  TH1F* rateHist_p[4];
  for(Int_t iter = 0; iter < 4; iter++){
    rateHist_p[iter] = new TH1F(Form("rateHist_%d_h", iter+2), Form("rateHist_%d_h", iter+2), 5, 30, 130);
  }


  HLTTree->SetBranchAddress("Event",&hlt_event);
  HLTTree->SetBranchAddress("Run",&hlt_run);
  HLTTree->SetBranchAddress("LumiBlock",&hlt_lumi);

  HLTTree->SetBranchAddress("L1IsolEmEt", L1IsolEmEt);
  HLTTree->SetBranchAddress("L1NIsolEmEt", L1NIsolEmEt);
  HLTTree->SetBranchAddress("L1CenJetEt", L1CenJetEt);
  HLTTree->SetBranchAddress("L1ForJetEt", L1ForJetEt);
  HLTTree->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7);
  HLTTree->SetBranchAddress("L1_SingleEG12", &L1_SingleEG12);
  HLTTree->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20);

  HLTTree->SetBranchAddress(ak2CaloName[0], &HLT_PuAK2CaloJet40_v1);
  HLTTree->SetBranchAddress(ak2CaloName[1], &HLT_PuAK2CaloJet60_v1);
  HLTTree->SetBranchAddress(ak2CaloName[2], &HLT_PuAK2CaloJet80_v1);
  HLTTree->SetBranchAddress(ak2CaloName[3], &HLT_PuAK2CaloJet100_v1);
  HLTTree->SetBranchAddress(ak2CaloName[4], &HLT_PuAK2CaloJet120_v1);

  HLTTree->SetBranchAddress(ak2CaloName_MID[0], &HLT_PuAK2CaloJet40_MID_v1);
  HLTTree->SetBranchAddress(ak2CaloName_MID[1], &HLT_PuAK2CaloJet60_MID_v1);
  HLTTree->SetBranchAddress(ak2CaloName_MID[2], &HLT_PuAK2CaloJet80_MID_v1);
  HLTTree->SetBranchAddress(ak2CaloName_MID[3], &HLT_PuAK2CaloJet100_MID_v1);
  HLTTree->SetBranchAddress(ak2CaloName_MID[4], &HLT_PuAK2CaloJet120_MID_v1);

  HLTTree->SetBranchAddress(ak3CaloName[0], &HLT_PuAK3CaloJet40_v1);
  HLTTree->SetBranchAddress(ak3CaloName[1], &HLT_PuAK3CaloJet60_v1);
  HLTTree->SetBranchAddress(ak3CaloName[2], &HLT_PuAK3CaloJet80_v1);
  HLTTree->SetBranchAddress(ak3CaloName[3], &HLT_PuAK3CaloJet100_v1);
  HLTTree->SetBranchAddress(ak3CaloName[4], &HLT_PuAK3CaloJet120_v1);

  HLTTree->SetBranchAddress(ak4CaloName[0], &HLT_PuAK4CaloJet40_v1);
  HLTTree->SetBranchAddress(ak4CaloName[1], &HLT_PuAK4CaloJet60_v1);
  HLTTree->SetBranchAddress(ak4CaloName[2], &HLT_PuAK4CaloJet80_v1);
  HLTTree->SetBranchAddress(ak4CaloName[3], &HLT_PuAK4CaloJet100_v1);
  HLTTree->SetBranchAddress(ak4CaloName[4], &HLT_PuAK4CaloJet120_v1);

  HLTTree->SetBranchAddress(ak4CaloName_MID[0], &HLT_PuAK4CaloJet40_MID_v1);
  HLTTree->SetBranchAddress(ak4CaloName_MID[1], &HLT_PuAK4CaloJet60_MID_v1);
  HLTTree->SetBranchAddress(ak4CaloName_MID[2], &HLT_PuAK4CaloJet80_MID_v1);
  HLTTree->SetBranchAddress(ak4CaloName_MID[3], &HLT_PuAK4CaloJet100_MID_v1);
  HLTTree->SetBranchAddress(ak4CaloName_MID[4], &HLT_PuAK4CaloJet120_MID_v1);

  HLTTree->SetBranchAddress(ak4CaloName_3Jet[0], &HLT_PuAK4CaloJet80_45_45_MID_v1);

  HLTTree->SetBranchAddress(ak4CaloName_ULTRAMID[0], &HLT_PuAK4CaloJet80_MID2_v1);
  HLTTree->SetBranchAddress(ak4CaloName_ULTRAMID[1], &HLT_PuAK4CaloJet80_2MID2_v1);
  HLTTree->SetBranchAddress(ak4CaloName_ULTRAMID[2], &HLT_PuAK4CaloJet80_MID3_v1);
  HLTTree->SetBranchAddress(ak4CaloName_ULTRAMID[3], &HLT_PuAK4CaloJet80_2MID3_v1);
  HLTTree->SetBranchAddress(ak4CaloName_ULTRAMID[4], &HLT_PuAK4CaloJet100_MID2_v1);
  HLTTree->SetBranchAddress(ak4CaloName_ULTRAMID[5], &HLT_PuAK4CaloJet100_2MID2_v1);
  HLTTree->SetBranchAddress(ak4CaloName_ULTRAMID[6], &HLT_PuAK4CaloJet100_MID3_v1);
  HLTTree->SetBranchAddress(ak4CaloName_ULTRAMID[7], &HLT_PuAK4CaloJet100_2MID3_v1);

  HLTTree->SetBranchAddress(ak5CaloName[0], &HLT_PuAK5CaloJet40_v1);
  HLTTree->SetBranchAddress(ak5CaloName[1], &HLT_PuAK5CaloJet60_v1);
  HLTTree->SetBranchAddress(ak5CaloName[2], &HLT_PuAK5CaloJet80_v1);
  HLTTree->SetBranchAddress(ak5CaloName[3], &HLT_PuAK5CaloJet100_v1);
  HLTTree->SetBranchAddress(ak5CaloName[4], &HLT_PuAK5CaloJet120_v1);

  HLTTree->SetBranchAddress(photonName[0], &HLT_HISinglePhoton10_v1);
  HLTTree->SetBranchAddress(photonName[1], &HLT_HISinglePhoton15_v1);
  HLTTree->SetBranchAddress(photonName[2], &HLT_HISinglePhoton20_v1);
  HLTTree->SetBranchAddress(photonName[3], &HLT_HISinglePhoton40_v1);
  HLTTree->SetBranchAddress(photonName[4], &HLT_HISinglePhoton60_v1);

  HLTTree->SetBranchAddress(photonName_MID[0], &HLT_HISinglePhoton10_MID_v1);
  HLTTree->SetBranchAddress(photonName_MID[1], &HLT_HISinglePhoton15_MID_v1);
  HLTTree->SetBranchAddress(photonName_MID[2], &HLT_HISinglePhoton20_MID_v1);
  HLTTree->SetBranchAddress(photonName_MID[3], &HLT_HISinglePhoton40_MID_v1);
  HLTTree->SetBranchAddress(photonName_MID[4], &HLT_HISinglePhoton60_MID_v1);

  HLTTree->SetBranchAddress(photonDoubleName[0], &HLT_HIDoublePhoton10_v1);
  HLTTree->SetBranchAddress(photonDoubleName[1], &HLT_HIDoublePhoton15_v1);
  HLTTree->SetBranchAddress(photonDoubleName[2], &HLT_HIDoublePhoton20_v1);
  HLTTree->SetBranchAddress(photonDoubleName[3], &HLT_HIDoublePhoton30_v1);
  HLTTree->SetBranchAddress(photonDoubleName[4], &HLT_HIPhoton15_Photon10_v1);
  HLTTree->SetBranchAddress(photonDoubleName[5], &HLT_HIPhoton20_Photon10_v1);
  HLTTree->SetBranchAddress(photonDoubleName[6], &HLT_HIPhoton20_Photon15_v1);


  TFile *AnaFile = TFile::Open(AnaFilename);
  TTree *Ana3CaloTree = (TTree*)AnaFile->Get(Ana3CaloTreename); 
  TTree *Ana4CaloTree = (TTree*)AnaFile->Get(Ana4CaloTreename); 
  TTree *AnaPhotonTree = (TTree*)AnaFile->Get(AnaPhotonTreename); 
  TTree *AnaHITree = (TTree*)AnaFile->Get(AnaHITreename); 
  TTree *AnaHLTTree = (TTree*)AnaFile->Get(AnaHLTTreename);
  TTree *AnaSkimTree = (TTree*)AnaFile->Get(AnaSkimTreename);
  TTree *AnaTrkTree = (TTree*)AnaFile->Get(AnaTrkTreename);
  TTree *AnaGenTree = (TTree*)AnaFile->Get(AnaGenTreename);

  Int_t ana_event, ana_lumi;//, ana_run, ana_lumi;
  Int_t hiBin;

  Int_t pcollisionEventSelection;

  Float_t L1IsolEmEt_ana[4];
  Float_t L1NIsolEmEt_ana[4];
  Float_t L1CenJetEt_ana[4];
  Float_t L1ForJetEt_ana[4];
  Int_t L1_SingleMu7_ana;
  Int_t L1_SingleEG12_ana;
  Int_t L1_SingleEG20_ana;

  Int_t nTrk;
  Bool_t trkFake[50000];
  Float_t trkPt[50000], trkPhi[50000], trkEta[50000];

  Int_t mult = 0;
  //  Int_t chg[100000];
  //  Float_t genPt[100000], genPhi[100000], genEta[100000];
  //  Float_t genPt[100000], genEta[100000];

  Int_t n3Caloref;
  Float_t jt3Calopt[500], jt3Caloeta[500], jt3Calophi[500];
  Float_t ref3Calopt[500], ref3Caloeta[500], ref3Calophi[500];
  Int_t n3Calogen;
  Float_t gen3Calopt[500], gen3Caloeta[500], gen3Calophi[500];

  Int_t n4Caloref;
  Float_t jt4Calopt[500], jt4Caloeta[500], jt4Calophi[500];
  Float_t ref4Calopt[500], ref4Caloeta[500], ref4Calophi[500];
  Int_t n4Calogen;
  Float_t gen4Calopt[500], gen4Caloeta[500], gen4Calophi[500];

  Int_t nPhotons;
  Float_t photonPt[50], photonEta[50], photonPhi[50];   //[nPhotons]

  AnaHITree->SetBranchStatus("*", 0);
  AnaHITree->SetBranchStatus("evt", 1);
  AnaHITree->SetBranchStatus("lumi", 1);
  AnaHITree->SetBranchStatus("hiBin", 1);

  AnaHITree->SetBranchAddress("evt", &ana_event);
  AnaHITree->SetBranchAddress("lumi", &ana_lumi);
  AnaHITree->SetBranchAddress("hiBin", &hiBin);

  AnaHLTTree->SetBranchStatus("*", 0);
  AnaHLTTree->SetBranchStatus("L1IsolEmEt", 1);
  AnaHLTTree->SetBranchStatus("L1NIsolEmEt", 1);
  AnaHLTTree->SetBranchStatus("L1CenJetEt", 1);
  AnaHLTTree->SetBranchStatus("L1ForJetEt", 1);
  AnaHLTTree->SetBranchStatus("L1_SingleMu7", 1);
  AnaHLTTree->SetBranchStatus("L1_SingleEG12", 1);
  AnaHLTTree->SetBranchStatus("L1_SingleEG20", 1);

  //  Ana3CaloTree->SetBranchAddress("run", &ana_run);
  AnaHLTTree->SetBranchAddress("L1IsolEmEt", L1IsolEmEt_ana);
  AnaHLTTree->SetBranchAddress("L1NIsolEmEt", L1NIsolEmEt_ana);
  AnaHLTTree->SetBranchAddress("L1CenJetEt", L1CenJetEt_ana);
  AnaHLTTree->SetBranchAddress("L1ForJetEt", L1ForJetEt_ana);
  AnaHLTTree->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7_ana);
  AnaHLTTree->SetBranchAddress("L1_SingleEG12", &L1_SingleEG12_ana);
  AnaHLTTree->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20_ana);


  AnaSkimTree->SetBranchStatus("*", 0);
  AnaSkimTree->SetBranchStatus("pcollisionEventSelection", 1);

  AnaSkimTree->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection);


  AnaTrkTree->SetBranchStatus("*", 0);
  AnaTrkTree->SetBranchStatus("nTrk", 1);
  AnaTrkTree->SetBranchStatus("trkFake", 1);
  AnaTrkTree->SetBranchStatus("trkPt", 1);
  AnaTrkTree->SetBranchStatus("trkPhi", 1);
  AnaTrkTree->SetBranchStatus("trkEta", 1);

  AnaTrkTree->SetBranchAddress("nTrk", &nTrk);
  AnaTrkTree->SetBranchAddress("trkFake", trkFake);
  AnaTrkTree->SetBranchAddress("trkPt", trkPt);
  AnaTrkTree->SetBranchAddress("trkPhi", trkPhi);
  AnaTrkTree->SetBranchAddress("trkEta", trkEta);

  AnaGenTree->SetBranchStatus("*", 0);
  
  AnaGenTree->SetBranchStatus("mult", 1);
  /*AnaGenTree->SetBranchStatus("chg", 1);
  AnaGenTree->SetBranchStatus("pt", 1);
  AnaGenTree->SetBranchStatus("phi", 1);
  AnaGenTree->SetBranchStatus("eta", 1);
  */
  AnaGenTree->SetBranchAddress("mult", &mult);
  /*  AnaGenTree->SetBranchAddress("chg", chg);
  AnaGenTree->SetBranchAddress("pt", genPt);
  AnaGenTree->SetBranchAddress("phi", genPhi);
  AnaGenTree->SetBranchAddress("eta", genEta);
  */

  Ana3CaloTree->SetBranchStatus("*", 0);
  Ana3CaloTree->SetBranchStatus("nref", 1);
  Ana3CaloTree->SetBranchStatus("jtpt", 1);
  Ana3CaloTree->SetBranchStatus("jteta", 1);
  Ana3CaloTree->SetBranchStatus("jtphi", 1);
  Ana3CaloTree->SetBranchStatus("refpt", 1);
  Ana3CaloTree->SetBranchStatus("refeta", 1);
  Ana3CaloTree->SetBranchStatus("refphi", 1);
  Ana3CaloTree->SetBranchStatus("ngen", 1);
  Ana3CaloTree->SetBranchStatus("genpt", 1);
  Ana3CaloTree->SetBranchStatus("genphi", 1);
  Ana3CaloTree->SetBranchStatus("geneta", 1);

  Ana3CaloTree->SetBranchAddress("nref", &n3Caloref);
  Ana3CaloTree->SetBranchAddress("jtpt", jt3Calopt);
  Ana3CaloTree->SetBranchAddress("jteta", jt3Caloeta);
  Ana3CaloTree->SetBranchAddress("jtphi", jt3Calophi);
  Ana3CaloTree->SetBranchAddress("refpt", ref3Calopt);
  Ana3CaloTree->SetBranchAddress("refeta", ref3Caloeta);
  Ana3CaloTree->SetBranchAddress("refphi", ref3Calophi);
  Ana3CaloTree->SetBranchAddress("ngen", &n3Calogen);
  Ana3CaloTree->SetBranchAddress("genpt", gen3Calopt);
  Ana3CaloTree->SetBranchAddress("genphi", gen3Calophi);
  Ana3CaloTree->SetBranchAddress("geneta", gen3Caloeta);

  Ana4CaloTree->SetBranchStatus("*", 0);
  Ana4CaloTree->SetBranchStatus("nref", 1);
  Ana4CaloTree->SetBranchStatus("jtpt", 1);
  Ana4CaloTree->SetBranchStatus("jteta", 1);
  Ana4CaloTree->SetBranchStatus("jtphi", 1);
  Ana4CaloTree->SetBranchStatus("refpt", 1);
  Ana4CaloTree->SetBranchStatus("refeta", 1);
  Ana4CaloTree->SetBranchStatus("refphi", 1);
  Ana4CaloTree->SetBranchStatus("ngen", 1);
  Ana4CaloTree->SetBranchStatus("genpt", 1);
  Ana4CaloTree->SetBranchStatus("genphi", 1);
  Ana4CaloTree->SetBranchStatus("geneta", 1);

  Ana4CaloTree->SetBranchAddress("nref", &n4Caloref);
  Ana4CaloTree->SetBranchAddress("jtpt", jt4Calopt);
  Ana4CaloTree->SetBranchAddress("jteta", jt4Caloeta);
  Ana4CaloTree->SetBranchAddress("jtphi", jt4Calophi);
  Ana4CaloTree->SetBranchAddress("refpt", ref4Calopt);
  Ana4CaloTree->SetBranchAddress("refeta", ref4Caloeta);
  Ana4CaloTree->SetBranchAddress("refphi", ref4Calophi);
  Ana4CaloTree->SetBranchAddress("ngen", &n4Calogen);
  Ana4CaloTree->SetBranchAddress("genpt", gen4Calopt);
  Ana4CaloTree->SetBranchAddress("genphi", gen4Calophi);
  Ana4CaloTree->SetBranchAddress("geneta", gen4Caloeta);

  AnaPhotonTree->SetBranchStatus("*", 0);
  AnaPhotonTree->SetBranchStatus("nPhotons", 1);
  AnaPhotonTree->SetBranchStatus("pt", 1);
  AnaPhotonTree->SetBranchStatus("eta", 1);
  AnaPhotonTree->SetBranchStatus("phi", 1);

  AnaPhotonTree->SetBranchAddress("nPhotons", &nPhotons);
  AnaPhotonTree->SetBranchAddress("pt", photonPt);
  AnaPhotonTree->SetBranchAddress("eta", photonEta);
  AnaPhotonTree->SetBranchAddress("phi", photonPhi);

  //book histos
  TH1D *hists2Calo_pt[n2Trig + 1], *hists2Calo_eta[n2Trig + 1];
  TH1D *hists3Calo_pt[n3Trig + 1], *hists3Calo_eta[n3Trig + 1];
  TH1D *hists4Calo_pt[n4Trig + 1], *hists4Calo_eta[n4Trig + 1];
  TH1D *hists5Calo_pt[n5Trig + 1], *hists5Calo_eta[n5Trig + 1];

  TH1D *hists2Calo_MID_pt[n2Trig + 1], *hists2Calo_MID_eta[n2Trig + 1];
  TH1D *hists4Calo_MID_pt[n4Trig + 1], *hists4Calo_MID_eta[n4Trig + 1];

  TH1D *hists4Calo_3Jet_pt[2], *hists4Calo_3Jet_eta[2];

  TH1D *hists2CaloTrk_pt[n2Trig + 1], *hists2CaloTrk_eta[n2Trig + 1];
  TH1D *hists3CaloTrk_pt[n3Trig + 1], *hists3CaloTrk_eta[n3Trig + 1];
  TH1D *hists4CaloTrk_pt[n4Trig + 1], *hists4CaloTrk_eta[n4Trig + 1];
  TH1D *hists5CaloTrk_pt[n5Trig + 1], *hists5CaloTrk_eta[n5Trig + 1];

  TH1D *hists2CaloTrk_MID_pt[n2Trig + 1], *hists2CaloTrk_MID_eta[n2Trig + 1];
  TH1D *hists4CaloTrk_MID_pt[n4Trig + 1], *hists4CaloTrk_MID_eta[n4Trig + 1];

  TH1D *hists2CaloGen_pt[n2Trig + 1], *hists2CaloGen_eta[n2Trig + 1];
  TH1D *hists3CaloGen_pt[n3Trig + 1], *hists3CaloGen_eta[n3Trig + 1];
  TH1D *hists4CaloGen_pt[n4Trig + 1], *hists4CaloGen_eta[n4Trig + 1];
  TH1D *hists5CaloGen_pt[n5Trig + 1], *hists5CaloGen_eta[n5Trig + 1];

  TH1D *histsPhoton_pt[nPhotonTrig + 1], *histsPhoton_eta[nPhotonTrig + 1];
  TH1D *histsPhoton_MID_pt[nPhotonTrig + 1], *histsPhoton_MID_eta[nPhotonTrig + 1];
  TH1D *histsDoublePhoton_pt[nPhotonTrig + 3], *histsDoublePhoton_eta[nPhotonTrig + 3];

  TH1I* hiBin_h = new TH1I("hiBin_h", "hiBin_h", 100, 0, 200);

  hists2Calo_pt[0] = new TH1D("leading_ak2Calo_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists2Calo_eta[0] = new TH1D("leading_ak2Calo_eta",";#eta^{jt}",nBins,-5,5);

  hists3Calo_pt[0] = new TH1D("leading_ak3Calo_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists3Calo_eta[0] = new TH1D("leading_ak3Calo_eta",";#eta^{jt}",nBins,-5,5);

  hists4Calo_pt[0] = new TH1D("leading_ak4Calo_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists4Calo_eta[0] = new TH1D("leading_ak4Calo_eta",";#eta^{jt}",nBins,-5,5);

  hists5Calo_pt[0] = new TH1D("leading_ak5Calo_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists5Calo_eta[0] = new TH1D("leading_ak5Calo_eta",";#eta^{jt}",nBins,-5,5);

  hists2Calo_MID_pt[0] = new TH1D("leading_ak2Calo_MID_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists2Calo_MID_eta[0] = new TH1D("leading_ak2Calo_MID_eta",";#eta^{jt}",nBins,-5,5);

  hists4Calo_MID_pt[0] = new TH1D("leading_ak4Calo_MID_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists4Calo_MID_eta[0] = new TH1D("leading_ak4Calo_MID_eta",";#eta^{jt}",nBins,-5,5);

  hists4Calo_3Jet_pt[0] = new TH1D("leading_ak4Calo_3Jet_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists4Calo_3Jet_eta[0] = new TH1D("leading_ak4Calo_3Jet_eta",";#eta^{jt}",nBins,-5,5);

  hists2CaloTrk_pt[0] = new TH1D("leading_ak2CaloTrk_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists2CaloTrk_eta[0] = new TH1D("leading_ak2CaloTrk_eta",";#eta^{jt}",nBins,-5,5);

  hists3CaloTrk_pt[0] = new TH1D("leading_ak3CaloTrk_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists3CaloTrk_eta[0] = new TH1D("leading_ak3CaloTrk_eta",";#eta^{jt}",nBins,-5,5);

  hists4CaloTrk_pt[0] = new TH1D("leading_ak4CaloTrk_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists4CaloTrk_eta[0] = new TH1D("leading_ak4CaloTrk_eta",";#eta^{jt}",nBins,-5,5);

  hists2CaloTrk_MID_pt[0] = new TH1D("leading_ak2CaloTrk_MID_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists2CaloTrk_MID_eta[0] = new TH1D("leading_ak2CaloTrk_MID_eta",";#eta^{jt}",nBins,-5,5);

  hists4CaloTrk_MID_pt[0] = new TH1D("leading_ak4CaloTrk_MID_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists4CaloTrk_MID_eta[0] = new TH1D("leading_ak4CaloTrk_MID_eta",";#eta^{jt}",nBins,-5,5);

  hists5CaloTrk_pt[0] = new TH1D("leading_ak5CaloTrk_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists5CaloTrk_eta[0] = new TH1D("leading_ak5CaloTrk_eta",";#eta^{jt}",nBins,-5,5);

  hists2CaloGen_pt[0] = new TH1D("leading_ak2CaloGen_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists2CaloGen_eta[0] = new TH1D("leading_ak2CaloGen_eta",";#eta^{jt}",nBins,-5,5);

  hists3CaloGen_pt[0] = new TH1D("leading_ak3CaloGen_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists3CaloGen_eta[0] = new TH1D("leading_ak3CaloGen_eta",";#eta^{jt}",nBins,-5,5);

  hists4CaloGen_pt[0] = new TH1D("leading_ak4CaloGen_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists4CaloGen_eta[0] = new TH1D("leading_ak4CaloGen_eta",";#eta^{jt}",nBins,-5,5);

  hists5CaloGen_pt[0] = new TH1D("leading_ak5CaloGen_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists5CaloGen_eta[0] = new TH1D("leading_ak5CaloGen_eta",";#eta^{jt}",nBins,-5,5);

  histsPhoton_pt[0] = new TH1D("leading_Photon_pt",";p_{T}^{jt}",nBins,0,maxpt);
  histsPhoton_eta[0] = new TH1D("leading_Photon_eta",";#eta^{jt}",nBins,-3,3);

  histsPhoton_MID_pt[0] = new TH1D("leading_Photon_MID_pt",";p_{T}^{jt}",nBins,0,maxpt);
  histsPhoton_MID_eta[0] = new TH1D("leading_Photon_MID_eta",";#eta^{jt}",nBins,-3,3);

  histsDoublePhoton_pt[0] = new TH1D("leading_DoublePhoton_pt",";p_{T}^{jt}",nBins,0,maxpt);
  histsDoublePhoton_eta[0] = new TH1D("leading_DoublePhoton_eta",";#eta^{jt}",nBins,-3,3);

  for(int i = 0; i < n2Trig; ++i){
    hists2Calo_pt[i+1] = (TH1D*)hists2Calo_pt[0]->Clone(ak2CaloName[i]);
    hists2Calo_eta[i+1] = (TH1D*)hists2Calo_eta[0]->Clone(ak2CaloName[i]+"eta");

    hists2Calo_MID_pt[i+1] = (TH1D*)hists2Calo_MID_pt[0]->Clone(ak2CaloName_MID[i]);
    hists2Calo_MID_eta[i+1] = (TH1D*)hists2Calo_MID_eta[0]->Clone(ak2CaloName_MID[i]+"eta");

    hists2CaloTrk_pt[i+1] = (TH1D*)hists2CaloTrk_pt[0]->Clone(ak2CaloTrkName[i]);
    hists2CaloTrk_eta[i+1] = (TH1D*)hists2CaloTrk_eta[0]->Clone(ak2CaloTrkName[i]+"eta");

    hists2CaloTrk_MID_pt[i+1] = (TH1D*)hists2CaloTrk_MID_pt[0]->Clone(ak2CaloTrkName_MID[i]);
    hists2CaloTrk_MID_eta[i+1] = (TH1D*)hists2CaloTrk_MID_eta[0]->Clone(ak2CaloTrkName_MID[i]+"eta");

    hists2CaloGen_pt[i+1] = (TH1D*)hists2CaloGen_pt[0]->Clone(ak2CaloGenName[i]);
    hists2CaloGen_eta[i+1] = (TH1D*)hists2CaloGen_eta[0]->Clone(ak2CaloGenName[i]+"eta");
  }

  for(int i = 0; i < n3Trig; ++i){
    hists3Calo_pt[i+1] = (TH1D*)hists3Calo_pt[0]->Clone(ak3CaloName[i]);
    hists3Calo_eta[i+1] = (TH1D*)hists3Calo_eta[0]->Clone(ak3CaloName[i]+"eta");

    hists3CaloTrk_pt[i+1] = (TH1D*)hists3CaloTrk_pt[0]->Clone(ak3CaloTrkName[i]);
    hists3CaloTrk_eta[i+1] = (TH1D*)hists3CaloTrk_eta[0]->Clone(ak3CaloTrkName[i]+"eta");

    hists3CaloGen_pt[i+1] = (TH1D*)hists3CaloGen_pt[0]->Clone(ak3CaloGenName[i]);
    hists3CaloGen_eta[i+1] = (TH1D*)hists3CaloGen_eta[0]->Clone(ak3CaloGenName[i]+"eta");
  }

  for(int i = 0; i < n4Trig; ++i){
    hists4Calo_pt[i+1] = (TH1D*)hists4Calo_pt[0]->Clone(ak4CaloName[i]);
    hists4Calo_eta[i+1] = (TH1D*)hists4Calo_eta[0]->Clone(ak4CaloName[i]+"eta");

    hists4Calo_MID_pt[i+1] = (TH1D*)hists4Calo_MID_pt[0]->Clone(ak4CaloName_MID[i]);
    hists4Calo_MID_eta[i+1] = (TH1D*)hists4Calo_MID_eta[0]->Clone(ak4CaloName_MID[i]+"eta");

    hists4CaloTrk_pt[i+1] = (TH1D*)hists4CaloTrk_pt[0]->Clone(ak4CaloTrkName[i]);
    hists4CaloTrk_eta[i+1] = (TH1D*)hists4CaloTrk_eta[0]->Clone(ak4CaloTrkName[i]+"eta");

    hists4CaloTrk_MID_pt[i+1] = (TH1D*)hists4CaloTrk_MID_pt[0]->Clone(ak4CaloTrkName_MID[i]);
    hists4CaloTrk_MID_eta[i+1] = (TH1D*)hists4CaloTrk_MID_eta[0]->Clone(ak4CaloTrkName_MID[i]+"eta");

    hists4CaloGen_pt[i+1] = (TH1D*)hists4CaloGen_pt[0]->Clone(ak4CaloGenName[i]);
    hists4CaloGen_eta[i+1] = (TH1D*)hists4CaloGen_eta[0]->Clone(ak4CaloGenName[i]+"eta");
  }

  hists4Calo_3Jet_pt[1] = (TH1D*)hists4Calo_3Jet_pt[0]->Clone(ak4CaloName_3Jet[0]);
  hists4Calo_3Jet_eta[1] = (TH1D*)hists4Calo_3Jet_eta[0]->Clone(ak4CaloName_3Jet[0]+"eta");

  for(int i = 0; i < n5Trig; ++i){
    hists5Calo_pt[i+1] = (TH1D*)hists5Calo_pt[0]->Clone(ak5CaloName[i]);
    hists5Calo_eta[i+1] = (TH1D*)hists5Calo_eta[0]->Clone(ak5CaloName[i]+"eta");

    hists5CaloTrk_pt[i+1] = (TH1D*)hists5CaloTrk_pt[0]->Clone(ak5CaloTrkName[i]);
    hists5CaloTrk_eta[i+1] = (TH1D*)hists5CaloTrk_eta[0]->Clone(ak5CaloTrkName[i]+"eta");

    hists5CaloGen_pt[i+1] = (TH1D*)hists5CaloGen_pt[0]->Clone(ak5CaloGenName[i]);
    hists5CaloGen_eta[i+1] = (TH1D*)hists5CaloGen_eta[0]->Clone(ak5CaloGenName[i]+"eta");
  }

  for(int i = 0; i < nPhotonTrig; ++i){
    histsPhoton_pt[i+1] = (TH1D*)histsPhoton_pt[0]->Clone(photonName[i]);
    histsPhoton_eta[i+1] = (TH1D*)histsPhoton_eta[0]->Clone(photonName[i]+"eta");

    histsPhoton_MID_pt[i+1] = (TH1D*)histsPhoton_MID_pt[0]->Clone(photonName_MID[i]);
    histsPhoton_MID_eta[i+1] = (TH1D*)histsPhoton_MID_eta[0]->Clone(photonName_MID[i]+"eta");
  }

  for(int i = 0; i < nPhotonTrig+2; ++i){
    histsDoublePhoton_pt[i+1] = (TH1D*)histsDoublePhoton_pt[0]->Clone(photonDoubleName[i]);
    histsDoublePhoton_eta[i+1] = (TH1D*)histsDoublePhoton_eta[0]->Clone(photonDoubleName[i]+"eta");
  }

  std::cout << "Events in HLT file: " << HLTTree->GetEntries() << std::endl;
  std::cout << "Events in Ana file: " << Ana3CaloTree->GetEntries() << std::endl;

  //make map
  EventMatchingCMS *matcher = new EventMatchingCMS();
  EventMatchingCMS *matcher2 = new EventMatchingCMS();

  std::vector<Int_t>* event_p = new std::vector<Int_t>;
  std::vector<Int_t>* jt_p = new std::vector<Int_t>;
  std::vector<Int_t>* trk_p = new std::vector<Int_t>;

  for(Long64_t entry = 0; entry < HLTTree->GetEntries(); ++entry)
  {
    HLTTree->GetEntry(entry);

    Int_t prematch3 = TMath::Abs(L1CenJetEt[0])/10;
    prematch3 *= 100;
    prematch3 += TMath::Abs(L1CenJetEt[1])/10;

    Int_t prematch3For = TMath::Abs(L1ForJetEt[0])/10;
    prematch3For *= 100;
    prematch3For += TMath::Abs(L1ForJetEt[1])/10;
    
    //    Int_t prematchEm = TMath::Abs(L1NIsolEmEt[0]);
    //    if(L1IsolEmEt[0] > 0) prematchEm += 1000;
    //    if(L1IsolEmEt[1] > 0) prematchEm += 100;
    
    //    Int_t match3 = prematch3;
    matcher->addEvent(hlt_event, hlt_lumi, 0, entry);
  }

  //  std::vector<Int_t>* caloEvt100_p = new std::vector<Int_t>;
  //  std::vector<Float_t>* caloPt100_p = new std::vector<Float_t>;

  std::cout << "6" << std::endl;

  // analysis loop
  std::cout << "COMMENCE LOOP" << std::endl;
  int matched = 0;
  for(Long64_t entry = 0; entry < Ana3CaloTree->GetEntries(); ++entry)
  {
    if(entry%10000 == 0) std::cout << entry << std::endl;

    Ana3CaloTree->GetEntry(entry);
    Ana4CaloTree->GetEntry(entry);
    AnaPhotonTree->GetEntry(entry);
    AnaHITree->GetEntry(entry);
    AnaHLTTree->GetEntry(entry);
    AnaSkimTree->GetEntry(entry);
    AnaTrkTree->GetEntry(entry);
    AnaGenTree->GetEntry(entry);

    Int_t prematch3_ana = TMath::Abs(L1CenJetEt_ana[0])/10;
    prematch3_ana *= 100;
    prematch3_ana += TMath::Abs(L1CenJetEt_ana[1])/10;

    Int_t prematch3For_ana = TMath::Abs(L1ForJetEt_ana[0])/10;
    prematch3For_ana *= 100;
    prematch3For_ana += TMath::Abs(L1ForJetEt_ana[1])/10;

    //    Int_t prematchEm_ana = TMath::Abs(L1NIsolEmEt_ana[0]);
    //    if(L1NIsolEmEt_ana[0] > 0) prematchEm_ana += 10;
    //    if(L1NIsolEmEt_ana[1] > 0) prematchEm_ana += 1;

    //    Int_t match3_ana = prematch3_ana;
    if(!pcollisionEventSelection) continue;

    long long hlt_entry = matcher->retrieveEvent(ana_event, ana_lumi, 0);
    if(hlt_entry == -1){
      matcher2->addEvent(ana_event, n4Caloref, nTrk, entry);
      continue;
    }
    HLTTree->GetEntry(hlt_entry);
    matched++;
    event_p->push_back(ana_event);
    jt_p->push_back(n4Caloref);
    trk_p->push_back(nTrk+mult);

    hiBin_h->Fill(hiBin);

    if(n3Caloref != 0){
      recoLeadJet3_p[0]->Fill(jt3Calopt[0]);
      recoLeadJet3_p[1]->Fill(jt3Caloeta[0]);
      recoLeadJet3_p[2]->Fill(jt3Calophi[0]);
    }

    if(n4Caloref != 0){
      recoLeadJet4_p[0]->Fill(jt4Calopt[0]);
      recoLeadJet4_p[1]->Fill(jt4Caloeta[0]);
      recoLeadJet4_p[2]->Fill(jt4Calophi[0]);
    }

    if(n3Calogen != 0){
      Float_t maxGenPt = -1;
      Float_t maxGenEta = -1;
      Float_t maxGenPhi = -1;
      
      for(Int_t iter = 0; iter < n3Calogen; iter++){
        if(gen3Calopt[iter] > maxGenPt){
          maxGenPt = gen3Calopt[iter];
          maxGenPhi = gen3Calophi[iter];
          maxGenEta = gen3Caloeta[iter];
        }
      }
      
      genLeadJet3_p[0]->Fill(maxGenPt);
      genLeadJet3_p[1]->Fill(maxGenEta);
      genLeadJet3_p[2]->Fill(maxGenPhi);
    }

    if(n4Calogen != 0){
      Float_t maxGenPt = -1;
      Float_t maxGenEta = -1;
      Float_t maxGenPhi = -1;
      
      for(Int_t iter = 0; iter < n4Calogen; iter++){
        if(gen4Calopt[iter] > maxGenPt){
          maxGenPt = gen4Calopt[iter];
          maxGenPhi = gen4Calophi[iter];
          maxGenEta = gen4Caloeta[iter];
        }
      }
      
      genLeadJet4_p[0]->Fill(maxGenPt);
      genLeadJet4_p[1]->Fill(maxGenEta);
      genLeadJet4_p[2]->Fill(maxGenPhi);
    }


    for(Int_t iter = 0; iter < n4Calogen; iter++){
      if(gen4Calopt[iter] > 40){
	nTrueJet40++;
	break;
      }
    }

    for(Int_t iter = 0; iter < n4Calogen; iter++){
      if(gen4Calopt[iter] > 80){
	nTrueJet80++;
	break;
      }
    }

    for(Int_t iter = 0; iter < n4Calogen; iter++){
      if(gen4Calopt[iter] > 90){
	nTrueJet90++;
	break;
      }
    }

    for(Int_t iter = 0; iter < n4Calogen; iter++){
      if(gen4Calopt[iter] > 100){
	nTrueJet100++;
	break;
      }
    }

    for(Int_t iter = 0; iter < n4Calogen; iter++){
      if(gen4Calopt[iter] > 110){
	nTrueJet110++;
	break;
      }
    }

    for(Int_t iter = 0; iter < n4Calogen; iter++){
      if(gen4Calopt[iter] > 120){
	nTrueJet120++;
	break;
      }
    }


    for(Int_t iter = 0; iter < n4Caloref; iter++){
      if(jt4Calopt[iter] > 40){
	nRecoJet40++;
	break;
      }
    }

    for(Int_t iter = 0; iter < n4Caloref; iter++){
      if(jt4Calopt[iter] > 80){
	nRecoJet80++;
	break;
      }
    }


    if(HLT_PuAK2CaloJet40_v1) n2CaloJet40++;
    if(HLT_PuAK2CaloJet60_v1) n2CaloJet60++;
    if(HLT_PuAK2CaloJet80_v1) n2CaloJet80++;
    if(HLT_PuAK2CaloJet100_v1) n2CaloJet100++;
    if(HLT_PuAK2CaloJet120_v1) n2CaloJet120++;

    if(HLT_PuAK2CaloJet40_MID_v1) n2CaloJet40_MID++;
    if(HLT_PuAK2CaloJet60_MID_v1) n2CaloJet60_MID++;
    if(HLT_PuAK2CaloJet80_MID_v1) n2CaloJet80_MID++;
    if(HLT_PuAK2CaloJet100_MID_v1) n2CaloJet100_MID++;
    if(HLT_PuAK2CaloJet120_MID_v1) n2CaloJet120_MID++;

    if(HLT_PuAK3CaloJet40_v1) n3CaloJet40++;
    if(HLT_PuAK3CaloJet60_v1) n3CaloJet60++;
    if(HLT_PuAK3CaloJet80_v1) n3CaloJet80++;
    if(HLT_PuAK3CaloJet100_v1) n3CaloJet100++;
    if(HLT_PuAK3CaloJet120_v1) n3CaloJet120++;

    if(HLT_PuAK4CaloJet40_v1) n4CaloJet40++;
    if(HLT_PuAK4CaloJet60_v1) n4CaloJet60++;
    if(HLT_PuAK4CaloJet80_v1) n4CaloJet80++;
    if(HLT_PuAK4CaloJet100_v1) n4CaloJet100++;
    if(HLT_PuAK4CaloJet120_v1) n4CaloJet120++;

    if(HLT_PuAK4CaloJet40_MID_v1) n4CaloJet40_MID++;
    if(HLT_PuAK4CaloJet60_MID_v1) n4CaloJet60_MID++;
    if(HLT_PuAK4CaloJet80_MID_v1) n4CaloJet80_MID++;
    if(HLT_PuAK4CaloJet100_MID_v1) n4CaloJet100_MID++;
    if(HLT_PuAK4CaloJet120_MID_v1) n4CaloJet120_MID++;

    if(HLT_PuAK4CaloJet80_45_45_MID_v1) n4CaloJet80_45_45_MID++;
    if(HLT_PuAK4CaloJet80_45_45_MID_v1 && !HLT_PuAK4CaloJet100_MID_v1) n4CaloJet80_45_45_MID_no100++;
    if(HLT_PuAK4CaloJet80_45_45_MID_v1 && !HLT_PuAK4CaloJet120_MID_v1) n4CaloJet80_45_45_MID_no120++;

    if(HLT_PuAK4CaloJet80_MID2_v1) n4CaloJet80_MID2++;
    if(HLT_PuAK4CaloJet80_2MID2_v1) n4CaloJet80_2MID2++;
    if(HLT_PuAK4CaloJet80_MID3_v1) n4CaloJet80_MID3++;
    if(HLT_PuAK4CaloJet80_2MID3_v1) n4CaloJet80_2MID3++;
    if(HLT_PuAK4CaloJet100_MID2_v1) n4CaloJet100_MID2++;
    if(HLT_PuAK4CaloJet100_2MID2_v1) n4CaloJet100_2MID2++;
    if(HLT_PuAK4CaloJet100_MID3_v1) n4CaloJet100_MID3++;
    if(HLT_PuAK4CaloJet100_2MID3_v1) n4CaloJet100_2MID3++;

    if(HLT_PuAK4CaloJet80_MID2_v1 && !HLT_PuAK4CaloJet100_MID_v1) n4CaloJet80_MID2_no100++;
    if(HLT_PuAK4CaloJet80_2MID2_v1 && !HLT_PuAK4CaloJet100_MID_v1) n4CaloJet80_2MID2_no100++;
    if(HLT_PuAK4CaloJet80_MID3_v1 && !HLT_PuAK4CaloJet100_MID_v1) n4CaloJet80_MID3_no100++;
    if(HLT_PuAK4CaloJet80_2MID3_v1 && !HLT_PuAK4CaloJet100_MID_v1) n4CaloJet80_2MID3_no100++;
    if(HLT_PuAK4CaloJet100_MID2_v1 && !HLT_PuAK4CaloJet120_MID_v1) n4CaloJet100_MID2_no120++;
    if(HLT_PuAK4CaloJet100_2MID2_v1 && !HLT_PuAK4CaloJet120_MID_v1) n4CaloJet100_2MID2_no120++;
    if(HLT_PuAK4CaloJet100_MID3_v1 && !HLT_PuAK4CaloJet120_MID_v1) n4CaloJet100_MID3_no120++;
    if(HLT_PuAK4CaloJet100_2MID3_v1 && !HLT_PuAK4CaloJet120_MID_v1) n4CaloJet100_2MID3_no120++;

    if(HLT_PuAK5CaloJet40_v1) n5CaloJet40++;
    if(HLT_PuAK5CaloJet60_v1) n5CaloJet60++;
    if(HLT_PuAK5CaloJet80_v1) n5CaloJet80++;
    if(HLT_PuAK5CaloJet100_v1) n5CaloJet100++;
    if(HLT_PuAK5CaloJet120_v1) n5CaloJet120++;

    if(HLT_PuAK4CaloJet100_v1 || HLT_HISinglePhoton40_v1) n4CaloJet100_SinglePhoton40++;
    if(HLT_PuAK4CaloJet100_MID_v1 || HLT_HISinglePhoton40_MID_v1) n4CaloJet100_SinglePhoton40_MID++;

    if(HLT_HISinglePhoton10_v1) nPhoton10++;
    if(HLT_HISinglePhoton15_v1) nPhoton15++;
    if(HLT_HISinglePhoton20_v1) nPhoton20++;
    if(HLT_HISinglePhoton40_v1) nPhoton40++;
    if(HLT_HISinglePhoton60_v1) nPhoton60++;

    if(HLT_HISinglePhoton10_MID_v1) nPhoton10_MID++;
    if(HLT_HISinglePhoton15_MID_v1) nPhoton15_MID++;
    if(HLT_HISinglePhoton20_MID_v1) nPhoton20_MID++;
    if(HLT_HISinglePhoton40_MID_v1) nPhoton40_MID++;
    if(HLT_HISinglePhoton60_MID_v1) nPhoton60_MID++;

    if(HLT_HIDoublePhoton10_v1) nDoublePhoton10++;
    if(HLT_HIDoublePhoton15_v1) nDoublePhoton15++;
    if(HLT_HIDoublePhoton20_v1) nDoublePhoton20++;
    if(HLT_HIDoublePhoton30_v1) nDoublePhoton30++;
    if(HLT_HIPhoton15_Photon10_v1) nPhoton15_10++;
    if(HLT_HIPhoton20_Photon10_v1) nPhoton20_10++;
    if(HLT_HIPhoton20_Photon15_v1) nPhoton20_15++;

    Double_t maxTrkPt = -1;
    Double_t maxTrkEta = -1;
  
    for(int i = 0; i < nTrk; i++){
      if(fabs(trkEta[i]) > 2.4) continue;
      if(trkFake[i]) continue;
      if(trkPt[i] > maxTrkPt){
	maxTrkPt = trkPt[i];
	maxTrkEta = trkEta[i];
      }
    }
  
    Double_t maxGenPt = -1;
    Double_t maxGenEta = -1;
    /*
    for(int i = 0; i < mult; i++){
      if(fabs(genEta[i]) > 2.4) continue;
      if(chg[i] == 0) continue;
      if(genPt[i] > maxGenPt){
	maxGenPt = genPt[i];
	maxGenEta = genEta[i];
      }
    }
    */
    Double_t max3CaloAnaPt = -1;
    Double_t max3CaloAnaEta = -100;
    
    for(int i = 0; i < n3Caloref; ++i){
      if(fabs(jt3Caloeta[i]) > 2.0) continue;
      if(jt3Calopt[i] > max3CaloAnaPt){
	max3CaloAnaPt = jt3Calopt[i];
	max3CaloAnaEta = jt3Caloeta[i];
      }
    }

    Double_t max4CaloAnaPt = -1;
    Double_t max4CaloAnaEta = -100;

    Double_t two4CaloAnaPt = -1;
    //    Double_t two4CaloAnaEta = -100;

    Double_t third4CaloAnaPt = -1;
    //    Double_t third4CaloAnaEta = -100;

    for(int i = 0; i < n4Caloref; ++i){
      if(fabs(jt4Caloeta[i]) > 2.0) continue;
      if(jt4Calopt[i] > max4CaloAnaPt){
	third4CaloAnaPt = two4CaloAnaPt;
	//        third4CaloAnaEta = two4CaloAnaEta;

	two4CaloAnaPt = max4CaloAnaPt;
	//	two4CaloAnaEta = max4CaloAnaEta;

	max4CaloAnaPt = jt4Calopt[i];
	max4CaloAnaEta = jt4Caloeta[i];
      }
      else if(jt4Calopt[i] > two4CaloAnaPt){
	third4CaloAnaPt = two4CaloAnaPt;
	//        third4CaloAnaEta = two4CaloAnaEta;

	two4CaloAnaPt = jt4Calopt[i];
	//	two4CaloAnaEta = jt4Caloeta[i];
      }
      else if(jt4Calopt[i] > third4CaloAnaPt){
	third4CaloAnaPt = jt4Calopt[i];
	//	third4CaloAnaEta = jt4Caloeta[i];
      }
    }

    Double_t max4GenAnaPt = -1;
    Double_t max4GenAnaEta = -100;

    for(int i = 0; i < n4Calogen; ++i){
      if(fabs(gen4Caloeta[i]) > 2.0) continue;
      if(gen4Calopt[i] > max4GenAnaPt){
	max4GenAnaPt = gen4Calopt[i];
	max4GenAnaEta = gen4Caloeta[i];
      }
    }

    if(max4CaloAnaPt > 100 && !HLT_PuAK4CaloJet100_v1){
      std::cout << "High Pt entry, jtpt, jteta, evt, lumi: " << entry << ", " << max4CaloAnaPt << ", " << max4CaloAnaEta << ", " << ana_event << ", " << ana_lumi << ", " << L1CenJetEt[0] << " ," << L1CenJetEt_ana[0] <<  std::endl;
    }

    if(max4CaloAnaPt > 80){
      totPt80_h->Fill(max4CaloAnaPt);
      totIntPt80_h->Fill(max4CaloAnaPt);
      if(max4CaloAnaPt > 85){
	totEta80_h->Fill(max4CaloAnaEta);
	totHiBin80_h->Fill(hiBin);
      }
    }
    if(max4CaloAnaPt > 100){
      totPt100_h->Fill(max4CaloAnaPt);
      totIntPt100_h->Fill(max4CaloAnaPt);
      if(max4CaloAnaPt > 105){
	totEta100_h->Fill(max4CaloAnaEta);
	totHiBin100_h->Fill(hiBin);
      }
    }
    if(max4CaloAnaPt > 120){
      totPt120_h->Fill(max4CaloAnaPt);
      totIntPt120_h->Fill(max4CaloAnaPt);
      if(max4CaloAnaPt > 125){
	totEta120_h->Fill(max4CaloAnaEta);
	totHiBin120_h->Fill(hiBin);
      }
    }

    Double_t maxPhotonAnaPt = -1;
    Double_t maxPhotonAnaEta = -100;

    Double_t twoPhotonAnaPt = -1;
    //    Double_t twoPhotonAnaEta = -100;
    for(int i = 0; i < nPhotons; ++i)
    {
      if(fabs(photonEta[i]) > 2.0) continue;
      if(photonPt[i] > maxPhotonAnaPt)
      {
	twoPhotonAnaPt = maxPhotonAnaPt;
	//	twoPhotonAnaEta = maxPhotonAnaEta;

	maxPhotonAnaPt = photonPt[i];
	maxPhotonAnaEta = photonEta[i];
      }
      else if(photonPt[i] > twoPhotonAnaPt){
        twoPhotonAnaPt = photonPt[i];
	//        twoPhotonAnaEta = photonEta[i];
      }
    }

    if(maxTrkPt > 0){
      hists2CaloTrk_pt[0]->Fill(maxTrkPt);
      hists2CaloTrk_eta[0]->Fill(maxTrkEta);

      hists2CaloTrk_MID_pt[0]->Fill(maxTrkPt);
      hists2CaloTrk_MID_eta[0]->Fill(maxTrkEta);

      hists3CaloTrk_pt[0]->Fill(maxTrkPt);
      hists3CaloTrk_eta[0]->Fill(maxTrkEta);

      hists4CaloTrk_pt[0]->Fill(maxTrkPt);
      hists4CaloTrk_eta[0]->Fill(maxTrkEta);

      hists4CaloTrk_MID_pt[0]->Fill(maxTrkPt);
      hists4CaloTrk_MID_eta[0]->Fill(maxTrkEta);

      hists5CaloTrk_pt[0]->Fill(maxTrkPt);
      hists5CaloTrk_eta[0]->Fill(maxTrkEta);

      if(HLT_PuAK2CaloJet40_v1){
        hists2CaloTrk_pt[1]->Fill(maxTrkPt);
        hists2CaloTrk_eta[1]->Fill(maxTrkEta);
      }
      if(HLT_PuAK2CaloJet60_v1){
        hists2CaloTrk_pt[2]->Fill(maxTrkPt);
        hists2CaloTrk_eta[2]->Fill(maxTrkEta);
      }
      if(HLT_PuAK2CaloJet80_v1){
        hists2CaloTrk_pt[3]->Fill(maxTrkPt);
        hists2CaloTrk_eta[3]->Fill(maxTrkEta);
      }
      if(HLT_PuAK2CaloJet100_v1){
        hists2CaloTrk_pt[4]->Fill(maxTrkPt);
        hists2CaloTrk_eta[4]->Fill(maxTrkEta);
      }
      if(HLT_PuAK2CaloJet120_v1){
        hists2CaloTrk_pt[5]->Fill(maxTrkPt);
        hists2CaloTrk_eta[5]->Fill(maxTrkEta);
      }

      if(HLT_PuAK2CaloJet40_MID_v1){
        hists2CaloTrk_MID_pt[1]->Fill(maxTrkPt);
        hists2CaloTrk_MID_eta[1]->Fill(maxTrkEta);
      }
      if(HLT_PuAK2CaloJet60_MID_v1){
        hists2CaloTrk_MID_pt[2]->Fill(maxTrkPt);
        hists2CaloTrk_MID_eta[2]->Fill(maxTrkEta);
      }
      if(HLT_PuAK2CaloJet80_MID_v1){
        hists2CaloTrk_MID_pt[3]->Fill(maxTrkPt);
        hists2CaloTrk_MID_eta[3]->Fill(maxTrkEta);
      }
      if(HLT_PuAK2CaloJet100_MID_v1){
        hists2CaloTrk_MID_pt[4]->Fill(maxTrkPt);
        hists2CaloTrk_MID_eta[4]->Fill(maxTrkEta);
      }
      if(HLT_PuAK2CaloJet120_MID_v1){
        hists2CaloTrk_MID_pt[5]->Fill(maxTrkPt);
        hists2CaloTrk_MID_eta[5]->Fill(maxTrkEta);
      }

      if(HLT_PuAK3CaloJet40_v1){
        hists3CaloTrk_pt[1]->Fill(maxTrkPt);
        hists3CaloTrk_eta[1]->Fill(maxTrkEta);
      }
      if(HLT_PuAK3CaloJet60_v1){
        hists3CaloTrk_pt[2]->Fill(maxTrkPt);
        hists3CaloTrk_eta[2]->Fill(maxTrkEta);
      }
      if(HLT_PuAK3CaloJet80_v1){
        hists3CaloTrk_pt[3]->Fill(maxTrkPt);
        hists3CaloTrk_eta[3]->Fill(maxTrkEta);
      }
      if(HLT_PuAK3CaloJet100_v1){
        hists3CaloTrk_pt[4]->Fill(maxTrkPt);
        hists3CaloTrk_eta[4]->Fill(maxTrkEta);
      }
      if(HLT_PuAK3CaloJet120_v1){
        hists3CaloTrk_pt[5]->Fill(maxTrkPt);
        hists3CaloTrk_eta[5]->Fill(maxTrkEta);
      }

      if(HLT_PuAK4CaloJet40_v1){
        hists4CaloTrk_pt[1]->Fill(maxTrkPt);
        hists4CaloTrk_eta[1]->Fill(maxTrkEta);
      }
      if(HLT_PuAK4CaloJet60_v1){
        hists4CaloTrk_pt[2]->Fill(maxTrkPt);
        hists4CaloTrk_eta[2]->Fill(maxTrkEta);
      }
      if(HLT_PuAK4CaloJet80_v1){
        hists4CaloTrk_pt[3]->Fill(maxTrkPt);
        hists4CaloTrk_eta[3]->Fill(maxTrkEta);
      }
      if(HLT_PuAK4CaloJet100_v1){
        hists4CaloTrk_pt[4]->Fill(maxTrkPt);
        hists4CaloTrk_eta[4]->Fill(maxTrkEta);
      }
      if(HLT_PuAK4CaloJet120_v1){
        hists4CaloTrk_pt[5]->Fill(maxTrkPt);
        hists4CaloTrk_eta[5]->Fill(maxTrkEta);
      }

      if(HLT_PuAK4CaloJet40_MID_v1){
        hists4CaloTrk_MID_pt[1]->Fill(maxTrkPt);
        hists4CaloTrk_MID_eta[1]->Fill(maxTrkEta);
      }
      if(HLT_PuAK4CaloJet60_MID_v1){
        hists4CaloTrk_MID_pt[2]->Fill(maxTrkPt);
        hists4CaloTrk_MID_eta[2]->Fill(maxTrkEta);
      }
      if(HLT_PuAK4CaloJet80_MID_v1){
        hists4CaloTrk_MID_pt[3]->Fill(maxTrkPt);
        hists4CaloTrk_MID_eta[3]->Fill(maxTrkEta);
      }
      if(HLT_PuAK4CaloJet100_MID_v1){
        hists4CaloTrk_MID_pt[4]->Fill(maxTrkPt);
        hists4CaloTrk_MID_eta[4]->Fill(maxTrkEta);
      }
      if(HLT_PuAK4CaloJet120_MID_v1){
        hists4CaloTrk_MID_pt[5]->Fill(maxTrkPt);
        hists4CaloTrk_MID_eta[5]->Fill(maxTrkEta);
      }

      if(HLT_PuAK5CaloJet40_v1){
        hists5CaloTrk_pt[1]->Fill(maxTrkPt);
        hists5CaloTrk_eta[1]->Fill(maxTrkEta);
      }
      if(HLT_PuAK5CaloJet60_v1){
        hists5CaloTrk_pt[2]->Fill(maxTrkPt);
        hists5CaloTrk_eta[2]->Fill(maxTrkEta);
      }
      if(HLT_PuAK5CaloJet80_v1){
        hists5CaloTrk_pt[3]->Fill(maxTrkPt);
        hists5CaloTrk_eta[3]->Fill(maxTrkEta);
      }
      if(HLT_PuAK5CaloJet100_v1){
        hists5CaloTrk_pt[4]->Fill(maxTrkPt);
        hists5CaloTrk_eta[4]->Fill(maxTrkEta);
      }
      if(HLT_PuAK5CaloJet120_v1){
        hists5CaloTrk_pt[5]->Fill(maxTrkPt);
        hists5CaloTrk_eta[5]->Fill(maxTrkEta);
      }
    }


    if(maxGenPt > 0){
      hists2CaloGen_pt[0]->Fill(maxGenPt);
      hists2CaloGen_eta[0]->Fill(maxGenEta);

      hists3CaloGen_pt[0]->Fill(maxGenPt);
      hists3CaloGen_eta[0]->Fill(maxGenEta);

      hists4CaloGen_pt[0]->Fill(maxGenPt);
      hists4CaloGen_eta[0]->Fill(maxGenEta);

      hists5CaloGen_pt[0]->Fill(maxGenPt);
      hists5CaloGen_eta[0]->Fill(maxGenEta);

      if(HLT_PuAK2CaloJet40_v1){
        hists2CaloGen_pt[1]->Fill(maxGenPt);
        hists2CaloGen_eta[1]->Fill(maxGenEta);
      }
      if(HLT_PuAK2CaloJet60_v1){
        hists2CaloGen_pt[2]->Fill(maxGenPt);
        hists2CaloGen_eta[2]->Fill(maxGenEta);
      }
      if(HLT_PuAK2CaloJet80_v1){
        hists2CaloGen_pt[3]->Fill(maxGenPt);
        hists2CaloGen_eta[3]->Fill(maxGenEta);
      }
      if(HLT_PuAK2CaloJet100_v1){
        hists2CaloGen_pt[4]->Fill(maxGenPt);
        hists2CaloGen_eta[4]->Fill(maxGenEta);
      }
      if(HLT_PuAK2CaloJet120_v1){
        hists2CaloGen_pt[5]->Fill(maxGenPt);
        hists2CaloGen_eta[5]->Fill(maxGenEta);
      }

      if(HLT_PuAK3CaloJet40_v1){
        hists3CaloGen_pt[1]->Fill(maxGenPt);
        hists3CaloGen_eta[1]->Fill(maxGenEta);
      }
      if(HLT_PuAK3CaloJet60_v1){
        hists3CaloGen_pt[2]->Fill(maxGenPt);
        hists3CaloGen_eta[2]->Fill(maxGenEta);
      }
      if(HLT_PuAK3CaloJet80_v1){
        hists3CaloGen_pt[3]->Fill(maxGenPt);
        hists3CaloGen_eta[3]->Fill(maxGenEta);
      }
      if(HLT_PuAK3CaloJet100_v1){
        hists3CaloGen_pt[4]->Fill(maxGenPt);
        hists3CaloGen_eta[4]->Fill(maxGenEta);
      }
      if(HLT_PuAK3CaloJet120_v1){
        hists3CaloGen_pt[5]->Fill(maxGenPt);
        hists3CaloGen_eta[5]->Fill(maxGenEta);
      }

      if(HLT_PuAK4CaloJet40_v1){
        hists4CaloGen_pt[1]->Fill(maxGenPt);
        hists4CaloGen_eta[1]->Fill(maxGenEta);
      }
      if(HLT_PuAK4CaloJet60_v1){
        hists4CaloGen_pt[2]->Fill(maxGenPt);
        hists4CaloGen_eta[2]->Fill(maxGenEta);
      }
      if(HLT_PuAK4CaloJet80_v1){
        hists4CaloGen_pt[3]->Fill(maxGenPt);
        hists4CaloGen_eta[3]->Fill(maxGenEta);
      }
      if(HLT_PuAK4CaloJet100_v1){
        hists4CaloGen_pt[4]->Fill(maxGenPt);
        hists4CaloGen_eta[4]->Fill(maxGenEta);
      }
      if(HLT_PuAK4CaloJet120_v1){
        hists4CaloGen_pt[5]->Fill(maxGenPt);
        hists4CaloGen_eta[5]->Fill(maxGenEta);
      }

      if(HLT_PuAK5CaloJet40_v1){
        hists5CaloGen_pt[1]->Fill(maxGenPt);
        hists5CaloGen_eta[1]->Fill(maxGenEta);
      }
      if(HLT_PuAK5CaloJet60_v1){
        hists5CaloGen_pt[2]->Fill(maxGenPt);
        hists5CaloGen_eta[2]->Fill(maxGenEta);
      }
      if(HLT_PuAK5CaloJet80_v1){
        hists5CaloGen_pt[3]->Fill(maxGenPt);
        hists5CaloGen_eta[3]->Fill(maxGenEta);
      }
      if(HLT_PuAK5CaloJet100_v1){
        hists5CaloGen_pt[4]->Fill(maxGenPt);
        hists5CaloGen_eta[4]->Fill(maxGenEta);
      }
      if(HLT_PuAK5CaloJet120_v1){
        hists5CaloGen_pt[5]->Fill(maxGenPt);
        hists5CaloGen_eta[5]->Fill(maxGenEta);
      }
    }

    if(max3CaloAnaPt > 0){
      hists2Calo_pt[0]->Fill(max3CaloAnaPt);
      hists2Calo_eta[0]->Fill(max3CaloAnaEta);

      hists2Calo_MID_pt[0]->Fill(max3CaloAnaPt);
      hists2Calo_MID_eta[0]->Fill(max3CaloAnaEta);

      hists3Calo_pt[0]->Fill(max3CaloAnaPt);
      hists3Calo_eta[0]->Fill(max3CaloAnaEta);

      if(HLT_PuAK2CaloJet40_v1){
	hists2Calo_pt[1]->Fill(max3CaloAnaPt);
	hists2Calo_eta[1]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK2CaloJet60_v1){
	hists2Calo_pt[2]->Fill(max3CaloAnaPt);
	hists2Calo_eta[2]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK2CaloJet80_v1){
	hists2Calo_pt[3]->Fill(max3CaloAnaPt);
	hists2Calo_eta[3]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK2CaloJet100_v1){
	hists2Calo_pt[4]->Fill(max3CaloAnaPt);
	hists2Calo_eta[4]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK2CaloJet120_v1){
	hists2Calo_pt[5]->Fill(max3CaloAnaPt);
	hists2Calo_eta[5]->Fill(max3CaloAnaEta);
      }

      if(HLT_PuAK2CaloJet40_MID_v1){
	hists2Calo_MID_pt[1]->Fill(max3CaloAnaPt);
	hists2Calo_MID_eta[1]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK2CaloJet60_MID_v1){
	hists2Calo_MID_pt[2]->Fill(max3CaloAnaPt);
	hists2Calo_MID_eta[2]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK2CaloJet80_MID_v1){
	hists2Calo_MID_pt[3]->Fill(max3CaloAnaPt);
	hists2Calo_MID_eta[3]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK2CaloJet100_MID_v1){
	hists2Calo_MID_pt[4]->Fill(max3CaloAnaPt);
	hists2Calo_MID_eta[4]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK2CaloJet120_MID_v1){
	hists2Calo_MID_pt[5]->Fill(max3CaloAnaPt);
	hists2Calo_MID_eta[5]->Fill(max3CaloAnaEta);
      }

      if(HLT_PuAK3CaloJet40_v1){
	hists3Calo_pt[1]->Fill(max3CaloAnaPt);
	hists3Calo_eta[1]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK3CaloJet60_v1){
	hists3Calo_pt[2]->Fill(max3CaloAnaPt);
	hists3Calo_eta[2]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK3CaloJet80_v1){
	hists3Calo_pt[3]->Fill(max3CaloAnaPt);
	hists3Calo_eta[3]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK3CaloJet100_v1){
	hists3Calo_pt[4]->Fill(max3CaloAnaPt);
	hists3Calo_eta[4]->Fill(max3CaloAnaEta);
      }
      if(HLT_PuAK3CaloJet120_v1){
	hists3Calo_pt[5]->Fill(max3CaloAnaPt);
	hists3Calo_eta[5]->Fill(max3CaloAnaEta);
      }
    }  
    else{
      if(HLT_PuAK2CaloJet40_v1) nFake2CaloJet40++;
      if(HLT_PuAK2CaloJet60_v1) nFake2CaloJet60++;
      if(HLT_PuAK2CaloJet80_v1) nFake2CaloJet80++;
      if(HLT_PuAK2CaloJet100_v1) nFake2CaloJet100++;
      if(HLT_PuAK2CaloJet120_v1) nFake2CaloJet120++;

      if(HLT_PuAK2CaloJet40_MID_v1) nFake2CaloJet40_MID++;
      if(HLT_PuAK2CaloJet60_MID_v1) nFake2CaloJet60_MID++;
      if(HLT_PuAK2CaloJet80_MID_v1) nFake2CaloJet80_MID++;
      if(HLT_PuAK2CaloJet100_MID_v1) nFake2CaloJet100_MID++;
      if(HLT_PuAK2CaloJet120_MID_v1) nFake2CaloJet120_MID++;

      if(HLT_PuAK3CaloJet40_v1) nFake3CaloJet40++;
      if(HLT_PuAK3CaloJet60_v1) nFake3CaloJet60++;
      if(HLT_PuAK3CaloJet80_v1) nFake3CaloJet80++;
      if(HLT_PuAK3CaloJet100_v1) nFake3CaloJet100++;
      if(HLT_PuAK3CaloJet120_v1) nFake3CaloJet120++;
    }

    if(max4CaloAnaPt > 0){
      hists4Calo_pt[0]->Fill(max4CaloAnaPt);
      hists4Calo_eta[0]->Fill(max4CaloAnaEta);

      hists4Calo_MID_pt[0]->Fill(max4CaloAnaPt);
      hists4Calo_MID_eta[0]->Fill(max4CaloAnaEta);

      if(two4CaloAnaPt > 55 && third4CaloAnaPt > 55){
	hists4Calo_3Jet_pt[0]->Fill(max4CaloAnaPt);
	hists4Calo_3Jet_eta[0]->Fill(max4CaloAnaEta);
      }

      hists5Calo_pt[0]->Fill(max4CaloAnaPt);
      hists5Calo_eta[0]->Fill(max4CaloAnaEta);

      if(HLT_PuAK4CaloJet40_v1){
	hists4Calo_pt[1]->Fill(max4CaloAnaPt);
	hists4Calo_eta[1]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK4CaloJet60_v1){
	hists4Calo_pt[2]->Fill(max4CaloAnaPt);
	hists4Calo_eta[2]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK4CaloJet80_v1){
	hists4Calo_pt[3]->Fill(max4CaloAnaPt);
	hists4Calo_eta[3]->Fill(max4CaloAnaEta);
      }
      else{
	if(max4CaloAnaPt > 80){
          maxMissedPt80_h->Fill(max4CaloAnaPt);
          maxMissedIntPt80_h->Fill(max4CaloAnaPt);
	  if(max4CaloAnaPt > 85){
	    maxMissedEta80_h->Fill(max4CaloAnaEta);
	    maxMissedHiBin80_h->Fill(hiBin);
	  }
	}
	if(max4CaloAnaPt > maxMissedPt80){
	  maxMissedPt80 = max4CaloAnaPt;
	  maxMissedEta80 = max4CaloAnaEta;
	}
	if(max4GenAnaPt > maxMissedGenPt80){
	  maxMissedGenPt80 = max4GenAnaPt;
	  maxMissedGenEta80 = max4GenAnaEta;
	}
      }
      if(HLT_PuAK4CaloJet100_v1){
	hists4Calo_pt[4]->Fill(max4CaloAnaPt);
	hists4Calo_eta[4]->Fill(max4CaloAnaEta);
      }
      else{
        if(max4CaloAnaPt > 100){
          maxMissedPt100_h->Fill(max4CaloAnaPt);
          maxMissedIntPt100_h->Fill(max4CaloAnaPt);
	  if(max4CaloAnaPt > 105){
	    maxMissedEta100_h->Fill(max4CaloAnaEta);
	    maxMissedHiBin100_h->Fill(hiBin);
	  }
	}
	if(max4CaloAnaPt > maxMissedPt100){
	  maxMissedPt100 = max4CaloAnaPt;
	  maxMissedEta100 = max4CaloAnaEta;
	}
	if(max4GenAnaPt > maxMissedGenPt100){
	  maxMissedGenPt100 = max4GenAnaPt;
	  maxMissedGenEta100 = max4GenAnaEta;
	}
      }
      if(HLT_PuAK4CaloJet120_v1){
	hists4Calo_pt[5]->Fill(max4CaloAnaPt);
	hists4Calo_eta[5]->Fill(max4CaloAnaEta);
      }
      else{
        if(max4CaloAnaPt > 120){
          maxMissedPt120_h->Fill(max4CaloAnaPt);
          maxMissedIntPt120_h->Fill(max4CaloAnaPt);
	  if(max4CaloAnaPt > 125){
	    maxMissedEta120_h->Fill(max4CaloAnaEta);
	    maxMissedHiBin120_h->Fill(hiBin);
	  }
	}
	if(max4CaloAnaPt > maxMissedPt120){
	  maxMissedPt120 = max4CaloAnaPt;
	  maxMissedEta120 = max4CaloAnaEta;
	}
	if(max4GenAnaPt > maxMissedGenPt120){
	  maxMissedGenPt120 = max4GenAnaPt;
	  maxMissedGenEta120 = max4GenAnaEta;
	}
      }

      if(HLT_PuAK4CaloJet40_MID_v1){
	hists4Calo_MID_pt[1]->Fill(max4CaloAnaPt);
	hists4Calo_MID_eta[1]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK4CaloJet60_MID_v1){
	hists4Calo_MID_pt[2]->Fill(max4CaloAnaPt);
	hists4Calo_MID_eta[2]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK4CaloJet80_MID_v1){
	hists4Calo_MID_pt[3]->Fill(max4CaloAnaPt);
	hists4Calo_MID_eta[3]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK4CaloJet100_MID_v1){
	hists4Calo_MID_pt[4]->Fill(max4CaloAnaPt);
	hists4Calo_MID_eta[4]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK4CaloJet120_MID_v1){
	hists4Calo_MID_pt[5]->Fill(max4CaloAnaPt);
	hists4Calo_MID_eta[5]->Fill(max4CaloAnaEta);
      }

      if(HLT_PuAK4CaloJet80_45_45_MID_v1 && two4CaloAnaPt > 55 && third4CaloAnaPt > 55){
	hists4Calo_3Jet_pt[1]->Fill(max4CaloAnaPt);
	hists4Calo_3Jet_eta[1]->Fill(max4CaloAnaEta);
      }

      if(HLT_PuAK5CaloJet40_v1){
	hists5Calo_pt[1]->Fill(max4CaloAnaPt);
	hists5Calo_eta[1]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK5CaloJet60_v1){
	hists5Calo_pt[2]->Fill(max4CaloAnaPt);
	hists5Calo_eta[2]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK5CaloJet80_v1){
	hists5Calo_pt[3]->Fill(max4CaloAnaPt);
	hists5Calo_eta[3]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK5CaloJet100_v1){
	hists5Calo_pt[4]->Fill(max4CaloAnaPt);
	hists5Calo_eta[4]->Fill(max4CaloAnaEta);
      }
      if(HLT_PuAK5CaloJet120_v1){
	hists5Calo_pt[5]->Fill(max4CaloAnaPt);
	hists5Calo_eta[5]->Fill(max4CaloAnaEta);
      }
    }  
    else{
      if(HLT_PuAK4CaloJet40_v1) nFake4CaloJet40++;
      if(HLT_PuAK4CaloJet60_v1) nFake4CaloJet60++;
      if(HLT_PuAK4CaloJet80_v1) nFake4CaloJet80++;
      if(HLT_PuAK4CaloJet100_v1) nFake4CaloJet100++;
      if(HLT_PuAK4CaloJet120_v1) nFake4CaloJet120++;

      if(HLT_PuAK4CaloJet40_MID_v1) nFake4CaloJet40_MID++;
      if(HLT_PuAK4CaloJet60_MID_v1) nFake4CaloJet60_MID++;
      if(HLT_PuAK4CaloJet80_MID_v1) nFake4CaloJet80_MID++;
      if(HLT_PuAK4CaloJet100_MID_v1) nFake4CaloJet100_MID++;
      if(HLT_PuAK4CaloJet120_MID_v1) nFake4CaloJet120_MID++;

      if(HLT_PuAK4CaloJet80_45_45_MID_v1) nFake4CaloJet80_45_45_MID++;
      if(HLT_PuAK4CaloJet80_45_45_MID_v1 && !HLT_PuAK4CaloJet100_MID_v1) nFake4CaloJet80_45_45_MID_no100++;
      if(HLT_PuAK4CaloJet80_45_45_MID_v1 && !HLT_PuAK4CaloJet120_MID_v1) nFake4CaloJet80_45_45_MID_no120++;

      if(HLT_PuAK5CaloJet40_v1) nFake5CaloJet40++;
      if(HLT_PuAK5CaloJet60_v1) nFake5CaloJet60++;
      if(HLT_PuAK5CaloJet80_v1) nFake5CaloJet80++;
      if(HLT_PuAK5CaloJet100_v1) nFake5CaloJet100++;
      if(HLT_PuAK5CaloJet120_v1) nFake5CaloJet120++;
    }
  


    if(maxPhotonAnaPt > 0){
      histsPhoton_pt[0]->Fill(maxPhotonAnaPt);
      histsPhoton_eta[0]->Fill(maxPhotonAnaEta);

      histsPhoton_MID_pt[0]->Fill(maxPhotonAnaPt);
      histsPhoton_MID_eta[0]->Fill(maxPhotonAnaEta);
      
      if(twoPhotonAnaPt > 5){
	histsDoublePhoton_pt[0]->Fill(maxPhotonAnaPt);
	histsDoublePhoton_eta[0]->Fill(maxPhotonAnaEta);
      }

      if(HLT_HISinglePhoton10_v1){
	histsPhoton_pt[1]->Fill(maxPhotonAnaPt);
	histsPhoton_eta[1]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HISinglePhoton15_v1){
	histsPhoton_pt[2]->Fill(maxPhotonAnaPt);
	histsPhoton_eta[2]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HISinglePhoton20_v1){
	histsPhoton_pt[3]->Fill(maxPhotonAnaPt);
	histsPhoton_eta[3]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HISinglePhoton40_v1){
	histsPhoton_pt[4]->Fill(maxPhotonAnaPt);
	histsPhoton_eta[4]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HISinglePhoton60_v1){
	histsPhoton_pt[5]->Fill(maxPhotonAnaPt);
	histsPhoton_eta[5]->Fill(maxPhotonAnaEta);
      }

      if(HLT_HISinglePhoton10_MID_v1){
	histsPhoton_MID_pt[1]->Fill(maxPhotonAnaPt);
	histsPhoton_MID_eta[1]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HISinglePhoton15_MID_v1){
	histsPhoton_MID_pt[2]->Fill(maxPhotonAnaPt);
	histsPhoton_MID_eta[2]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HISinglePhoton20_MID_v1){
	histsPhoton_MID_pt[3]->Fill(maxPhotonAnaPt);
	histsPhoton_MID_eta[3]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HISinglePhoton40_MID_v1){
	histsPhoton_MID_pt[4]->Fill(maxPhotonAnaPt);
	histsPhoton_MID_eta[4]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HISinglePhoton60_MID_v1){
	histsPhoton_MID_pt[5]->Fill(maxPhotonAnaPt);
	histsPhoton_MID_eta[5]->Fill(maxPhotonAnaEta);
      }

      if(HLT_HIDoublePhoton10_v1 && twoPhotonAnaPt > 5){
	histsDoublePhoton_pt[1]->Fill(maxPhotonAnaPt);
	histsDoublePhoton_eta[1]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HIDoublePhoton15_v1 && twoPhotonAnaPt > 0){
	histsDoublePhoton_pt[2]->Fill(maxPhotonAnaPt);
	histsDoublePhoton_eta[2]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HIDoublePhoton20_v1 && twoPhotonAnaPt > 0){
	histsDoublePhoton_pt[3]->Fill(maxPhotonAnaPt);
	histsDoublePhoton_eta[3]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HIDoublePhoton30_v1 && twoPhotonAnaPt > 0){
	histsDoublePhoton_pt[4]->Fill(maxPhotonAnaPt);
	histsDoublePhoton_eta[4]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HIPhoton15_Photon10_v1 && twoPhotonAnaPt > 0){
	histsDoublePhoton_pt[5]->Fill(maxPhotonAnaPt);
	histsDoublePhoton_eta[5]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HIPhoton20_Photon10_v1 && twoPhotonAnaPt > 0){
	histsDoublePhoton_pt[6]->Fill(maxPhotonAnaPt);
	histsDoublePhoton_eta[6]->Fill(maxPhotonAnaEta);
      }
      if(HLT_HIPhoton20_Photon15_v1 && twoPhotonAnaPt > 0){
	histsDoublePhoton_pt[7]->Fill(maxPhotonAnaPt);
	histsDoublePhoton_eta[7]->Fill(maxPhotonAnaEta);
      }
    }  
    else{
      if(HLT_HISinglePhoton10_v1) nFakePhoton10++;
      if(HLT_HISinglePhoton15_v1) nFakePhoton15++;
      if(HLT_HISinglePhoton20_v1) nFakePhoton20++;
      if(HLT_HISinglePhoton40_v1) nFakePhoton40++;
      if(HLT_HISinglePhoton60_v1){
	nFakePhoton60++;
	std::cout << ana_event << std::endl;
      }


      if(HLT_HISinglePhoton10_MID_v1) nFakePhoton10_MID++;
      if(HLT_HISinglePhoton15_MID_v1) nFakePhoton15_MID++;
      if(HLT_HISinglePhoton20_MID_v1) nFakePhoton20_MID++;
      if(HLT_HISinglePhoton40_MID_v1) nFakePhoton40_MID++;
      if(HLT_HISinglePhoton60_MID_v1) nFakePhoton60_MID++;

      if(HLT_HIDoublePhoton10_v1) nFakeDoublePhoton10++;
      if(HLT_HIDoublePhoton15_v1) nFakeDoublePhoton15++;
      if(HLT_HIDoublePhoton20_v1) nFakeDoublePhoton20++;
      if(HLT_HIDoublePhoton30_v1) nFakeDoublePhoton30++;
      if(HLT_HIPhoton15_Photon10_v1) nFakePhoton15_10++;
      if(HLT_HIPhoton20_Photon10_v1) nFakePhoton20_10++;
      if(HLT_HIPhoton20_Photon15_v1) nFakePhoton20_15++;
    }
  }


  std::cout << "Events matched: " << matched << std::endl;

  std::cout << "True jet > 40: " << nTrueJet40 << std::endl;
  std::cout << "True jet > 80: " << nTrueJet80 << std::endl;
  std::cout << "True jet > 90: " << nTrueJet90 << std::endl;
  std::cout << "True jet > 100: " << nTrueJet100 << std::endl;
  std::cout << "True jet > 110: " << nTrueJet110 << std::endl;
  std::cout << "True jet > 120: " << nTrueJet120 << std::endl;


  std::cout << "Reco jet > 40: " << nRecoJet40 << std::endl;
  std::cout << "Reco jet > 80: " << nRecoJet80 << std::endl;

  std::cout << "True trigger rates: " << std::endl;

  std::cout << "HLT_PuAK2CaloJet40_v1: " << n2CaloJet40 << ", " << nFake2CaloJet40 << std::endl;
  std::cout << "HLT_PuAK2CaloJet60_v1: " << n2CaloJet60 << ", " << nFake2CaloJet60 << std::endl;
  std::cout << "HLT_PuAK2CaloJet80_v1: " << n2CaloJet80 << ", " << nFake2CaloJet80 << std::endl;
  std::cout << "HLT_PuAK2CaloJet100_v1: " << n2CaloJet100 << ", " << nFake2CaloJet100 << std::endl;
  std::cout << "HLT_PuAK2CaloJet120_v1: " << n2CaloJet120 << ", " << nFake2CaloJet120 << std::endl;

  rateHist_p[0]->SetBinContent(1, n2CaloJet40);
  rateHist_p[0]->SetBinContent(2, n2CaloJet60);
  rateHist_p[0]->SetBinContent(3, n2CaloJet80);
  rateHist_p[0]->SetBinContent(4, n2CaloJet100);
  rateHist_p[0]->SetBinContent(5, n2CaloJet120);

  rateHist_p[0]->SetBinError(1, TMath::Sqrt(n2CaloJet40));
  rateHist_p[0]->SetBinError(2, TMath::Sqrt(n2CaloJet60));
  rateHist_p[0]->SetBinError(3, TMath::Sqrt(n2CaloJet80));
  rateHist_p[0]->SetBinError(4, TMath::Sqrt(n2CaloJet100));
  rateHist_p[0]->SetBinError(5, TMath::Sqrt(n2CaloJet120));

  for(Int_t iter = 0; iter < 5; iter++){
    rateHist_p[0]->SetBinContent(iter+1, rateHist_p[0]->GetBinContent(iter+1)*30000./matched);
    rateHist_p[0]->SetBinError(iter+1, rateHist_p[0]->GetBinError(iter+1)*30000./matched);
  }


  std::cout << "HLT_PuAK3CaloJet40_v1: " << n3CaloJet40 << ", " << nFake3CaloJet40 << std::endl;
  std::cout << "HLT_PuAK3CaloJet60_v1: " << n3CaloJet60 << ", " << nFake3CaloJet60 << std::endl;
  std::cout << "HLT_PuAK3CaloJet80_v1: " << n3CaloJet80 << ", " << nFake3CaloJet80 << std::endl;
  std::cout << "HLT_PuAK3CaloJet100_v1: " << n3CaloJet100 << ", " << nFake3CaloJet100 << std::endl;
  std::cout << "HLT_PuAK3CaloJet120_v1: " << n3CaloJet120 << ", " << nFake3CaloJet120 << std::endl;

  rateHist_p[1]->SetBinContent(1, n3CaloJet40);
  rateHist_p[1]->SetBinContent(2, n3CaloJet60);
  rateHist_p[1]->SetBinContent(3, n3CaloJet80);
  rateHist_p[1]->SetBinContent(4, n3CaloJet100);
  rateHist_p[1]->SetBinContent(5, n3CaloJet120);

  rateHist_p[1]->SetBinError(1, TMath::Sqrt(n3CaloJet40));
  rateHist_p[1]->SetBinError(2, TMath::Sqrt(n3CaloJet60));
  rateHist_p[1]->SetBinError(3, TMath::Sqrt(n3CaloJet80));
  rateHist_p[1]->SetBinError(4, TMath::Sqrt(n3CaloJet100));
  rateHist_p[1]->SetBinError(5, TMath::Sqrt(n3CaloJet120));

  for(Int_t iter = 0; iter < 5; iter++){
    rateHist_p[1]->SetBinContent(iter+1, rateHist_p[1]->GetBinContent(iter+1)*30000./matched);
    rateHist_p[1]->SetBinError(iter+1, rateHist_p[1]->GetBinError(iter+1)*30000./matched);
  }

  std::cout << "HLT_PuAK4CaloJet40_v1: " << n4CaloJet40 << ", " << nFake4CaloJet40 << std::endl;
  std::cout << "HLT_PuAK4CaloJet60_v1: " << n4CaloJet60 << ", " << nFake4CaloJet60 << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_v1: " << n4CaloJet80 << ", " << nFake4CaloJet80 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_v1: " << n4CaloJet100 << ", " << nFake4CaloJet100 << std::endl;
  std::cout << "HLT_PuAK4CaloJet120_v1: " << n4CaloJet120 << ", " << nFake4CaloJet120 << std::endl;

  std::cout << "Max missed pt 80, 100, 120: " << maxMissedPt80 << ", " << maxMissedPt100 << ", " << maxMissedPt120 << std::endl;
  std::cout << "Max missed eta 80, 100, 120: " << maxMissedEta80 << ", " << maxMissedEta100 << ", " << maxMissedEta120 << std::endl;

  std::cout << "Max missed gen pt 80, 100, 120: " << maxMissedGenPt80 << ", " << maxMissedGenPt100 << ", " << maxMissedGenPt120 << std::endl;
  std::cout << "Max missed gen eta 80, 100, 120: " << maxMissedGenEta80 << ", " << maxMissedGenEta100 << ", " << maxMissedGenEta120 << std::endl;

  std::cout << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_v1 && photon40: " << n4CaloJet100_SinglePhoton40 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_v1 && photon40 MID: " << n4CaloJet100_SinglePhoton40_MID << std::endl;
  std::cout << std::endl;

  rateHist_p[2]->SetBinContent(1, n4CaloJet40);
  rateHist_p[2]->SetBinContent(2, n4CaloJet60);
  rateHist_p[2]->SetBinContent(3, n4CaloJet80);
  rateHist_p[2]->SetBinContent(4, n4CaloJet100);
  rateHist_p[2]->SetBinContent(5, n4CaloJet120);

  rateHist_p[2]->SetBinError(1, TMath::Sqrt(n4CaloJet40));
  rateHist_p[2]->SetBinError(2, TMath::Sqrt(n4CaloJet60));
  rateHist_p[2]->SetBinError(3, TMath::Sqrt(n4CaloJet80));
  rateHist_p[2]->SetBinError(4, TMath::Sqrt(n4CaloJet100));
  rateHist_p[2]->SetBinError(5, TMath::Sqrt(n4CaloJet120));

  for(Int_t iter = 0; iter < 5; iter++){
    rateHist_p[2]->SetBinContent(iter+1, rateHist_p[2]->GetBinContent(iter+1)*30000./matched);
    rateHist_p[2]->SetBinError(iter+1, rateHist_p[2]->GetBinError(iter+1)*30000./matched);
  }


  std::cout << "HLT_PuAK5CaloJet40_v1: " << n5CaloJet40 << ", " << nFake5CaloJet40 << std::endl;
  std::cout << "HLT_PuAK5CaloJet60_v1: " << n5CaloJet60 << ", " << nFake5CaloJet60 << std::endl;
  std::cout << "HLT_PuAK5CaloJet80_v1: " << n5CaloJet80 << ", " << nFake5CaloJet80 << std::endl;
  std::cout << "HLT_PuAK5CaloJet100_v1: " << n5CaloJet100 << ", " << nFake5CaloJet100 << std::endl;
  std::cout << "HLT_PuAK5CaloJet120_v1: " << n5CaloJet120 << ", " << nFake5CaloJet120 << std::endl;

  rateHist_p[3]->SetBinContent(1, n5CaloJet40);
  rateHist_p[3]->SetBinContent(2, n5CaloJet60);
  rateHist_p[3]->SetBinContent(3, n5CaloJet80);
  rateHist_p[3]->SetBinContent(4, n5CaloJet100);
  rateHist_p[3]->SetBinContent(5, n5CaloJet120);

  rateHist_p[3]->SetBinError(1, TMath::Sqrt(n5CaloJet40));
  rateHist_p[3]->SetBinError(2, TMath::Sqrt(n5CaloJet60));
  rateHist_p[3]->SetBinError(3, TMath::Sqrt(n5CaloJet80));
  rateHist_p[3]->SetBinError(4, TMath::Sqrt(n5CaloJet100));
  rateHist_p[3]->SetBinError(5, TMath::Sqrt(n5CaloJet120));

  for(Int_t iter = 0; iter < 5; iter++){
    rateHist_p[3]->SetBinContent(iter+1, rateHist_p[3]->GetBinContent(iter+1)*30000./matched);
    rateHist_p[3]->SetBinError(iter+1, rateHist_p[3]->GetBinError(iter+1)*30000./matched);
  }

  std::cout << "HLT_PuAK2CaloJet40_MID_v1: " << n2CaloJet40_MID << ", " << nFake2CaloJet40_MID << std::endl;
  std::cout << "HLT_PuAK2CaloJet60_MID_v1: " << n2CaloJet60_MID << ", " << nFake2CaloJet60_MID << std::endl;
  std::cout << "HLT_PuAK2CaloJet80_MID_v1: " << n2CaloJet80_MID << ", " << nFake2CaloJet80_MID << std::endl;
  std::cout << "HLT_PuAK2CaloJet100_MID_v1: " << n2CaloJet100_MID << ", " << nFake2CaloJet100_MID << std::endl;
  std::cout << "HLT_PuAK2CaloJet120_MID_v1: " << n2CaloJet120_MID << ", " << nFake2CaloJet120_MID << std::endl;


  std::cout << "HLT_PuAK4CaloJet40_MID_v1: " << n4CaloJet40_MID << ", " << nFake4CaloJet40_MID << std::endl;
  std::cout << "HLT_PuAK4CaloJet60_MID_v1: " << n4CaloJet60_MID << ", " << nFake4CaloJet60_MID << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_MID_v1: " << n4CaloJet80_MID << ", " << nFake4CaloJet80_MID << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_MID_v1: " << n4CaloJet100_MID << ", " << nFake4CaloJet100_MID << std::endl;
  std::cout << "HLT_PuAK4CaloJet120_MID_v1: " << n4CaloJet120_MID << ", " << nFake4CaloJet120_MID << std::endl;

  std::cout << "HLT_PuAK4CaloJet80_45_45_v1: " << n4CaloJet80_45_45_MID << ", " << nFake4CaloJet80_45_45_MID << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_45_45_v1, no 100: " << n4CaloJet80_45_45_MID_no100 << ", " << nFake4CaloJet80_45_45_MID_no100 << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_45_45_v1, no 120: " << n4CaloJet80_45_45_MID_no120 << ", " << nFake4CaloJet80_45_45_MID_no120 << std::endl;

  std::cout << "HLT_PuAK4CaloJet80_MID2_v1: " << n4CaloJet80_MID2 << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_2MID2_v1: " << n4CaloJet80_2MID2 << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_MID3_v1: " << n4CaloJet80_MID3 << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_2MID3_v1: " << n4CaloJet80_2MID3 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_MID2_v1: " << n4CaloJet100_MID2 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_2MID2_v1: " << n4CaloJet100_2MID2 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_MID3_v1: " << n4CaloJet100_MID3 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_2MID3_v1: " << n4CaloJet100_2MID3 << std::endl;



  std::cout << "HLT_PuAK4CaloJet80_MID2_no100_v1: " << n4CaloJet80_MID2_no100 << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_2MID2_no100_v1: " << n4CaloJet80_2MID2_no100 << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_MID3_no100_v1: " << n4CaloJet80_MID3_no100 << std::endl;
  std::cout << "HLT_PuAK4CaloJet80_2MID3_no100_v1: " << n4CaloJet80_2MID3_no100 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_MID2_no120_v1: " << n4CaloJet100_MID2_no120 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_2MID2_no120_v1: " << n4CaloJet100_2MID2_no120 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_MID3_no120_v1: " << n4CaloJet100_MID3_no120 << std::endl;
  std::cout << "HLT_PuAK4CaloJet100_2MID3_no120_v1: " << n4CaloJet100_2MID3_no120 << std::endl;

  std::cout << "HLT_HISinglePhoton10_v1: " << nPhoton10 << ", " << nFakePhoton10 << std::endl;
  std::cout << "HLT_HISinglePhoton15_v1: " << nPhoton15 << ", " << nFakePhoton15 << std::endl;
  std::cout << "HLT_HISinglePhoton20_v1: " << nPhoton20 << ", " << nFakePhoton20 << std::endl;
  std::cout << "HLT_HISinglePhoton40_v1: " << nPhoton40 << ", " << nFakePhoton40 << std::endl;
  std::cout << "HLT_HISinglePhoton60_v1: " << nPhoton60 << ", " << nFakePhoton60 << std::endl;

  std::cout << "HLT_HISinglePhoton10_MID_v1: " << nPhoton10_MID << ", " << nFakePhoton10_MID << std::endl;
  std::cout << "HLT_HISinglePhoton15_MID_v1: " << nPhoton15_MID << ", " << nFakePhoton15_MID << std::endl;
  std::cout << "HLT_HISinglePhoton20_MID_v1: " << nPhoton20_MID << ", " << nFakePhoton20_MID << std::endl;
  std::cout << "HLT_HISinglePhoton40_MID_v1: " << nPhoton40_MID << ", " << nFakePhoton40_MID << std::endl;
  std::cout << "HLT_HISinglePhoton60_MID_v1: " << nPhoton60_MID << ", " << nFakePhoton60_MID << std::endl;

  std::cout << "HLT_HIDoublePhoton10_v1: " << nDoublePhoton10 << ", " << nFakeDoublePhoton10 << std::endl;
  std::cout << "HLT_HIDoublePhoton15_v1: " << nDoublePhoton15 << ", " << nFakeDoublePhoton15 << std::endl;
  std::cout << "HLT_HIDoublePhoton20_v1: " << nDoublePhoton20 << ", " << nFakeDoublePhoton20 << std::endl;
  std::cout << "HLT_HIDoublePhoton30_v1: " << nDoublePhoton30 << ", " << nFakeDoublePhoton30 << std::endl;

  std::cout << "HLT_HIPhoton15_Photon10_v1: " << nPhoton15_10 << ", " << nFakePhoton15_10 << std::endl;
  std::cout << "HLT_HIPhoton20_Photon10_v1: " << nPhoton20_10 << ", " << nFakePhoton20_10 << std::endl;
  std::cout << "HLT_HIPhoton20_Photon15_v1: " << nPhoton20_15 << ", " << nFakePhoton20_15 << std::endl;

  std::cout << std::endl;

  Int_t matched2 = 0;
  for(Int_t iter = 0; iter < matched; iter++){
    Int_t test = matcher2->retrieveEvent(event_p->at(iter), jt_p->at(iter), trk_p->at(iter));
    if(test == -1) continue;
    matched2++;
    std::cout << "event: " << event_p->at(iter) << std::endl;
  }

  std::cout <<"matched2: " <<matched2 << std::endl;


  /*
  std::cout << "Check 100 triggers and pt" << std::endl;

  std::cout << std::endl;
  Int_t nMatch_ = 0;
  for(Int_t pfIter = 0; pfIter < (Int_t)pfEvt100_p->size(); pfIter++){
    for(Int_t caloIter = 0; caloIter < (Int_t)caloEvt100_p->size(); caloIter++){

      if(caloEvt100_p->at(caloIter) == pfEvt100_p->at(pfIter)){
	std::cout << "Match #, entry #, caloPt, pfPt: " << nMatch_ << ", " << caloEvt100_p->at(caloIter) << ", " << caloPt100_p->at(caloIter) << ", " << pfPt100_p->at(pfIter) << std::endl;

	if(caloPt100_p->at(caloIter) - pfPt100_p->at(pfIter) > 0) std::cout << "SPECIAL SNOWFLAKE" << std::endl;
	nMatch_++;
	break;
      }
    }
  }

  std::cout << std::endl;
  std::cout << "nMatch: " << nMatch_ << std::endl;
  */
 
  //make turn-on curves
  TGraphAsymmErrors *a2Calo_pt[n2Trig], *a2Calo_eta[n2Trig];
  TGraphAsymmErrors *a3Calo_pt[n3Trig], *a3Calo_eta[n3Trig];
  TGraphAsymmErrors *a4Calo_pt[n4Trig], *a4Calo_eta[n4Trig];
  TGraphAsymmErrors *a5Calo_pt[n5Trig], *a5Calo_eta[n5Trig];

  TGraphAsymmErrors *a2Calo_MID_pt[n2Trig], *a2Calo_MID_eta[n2Trig];
  TGraphAsymmErrors *a4Calo_MID_pt[n4Trig], *a4Calo_MID_eta[n4Trig];
  TGraphAsymmErrors *a4Calo_3Jet_pt[n4Trig], *a4Calo_3Jet_eta[n4Trig];

  TGraphAsymmErrors *a2CaloTrk_pt[n2Trig], *a2CaloTrk_eta[n2Trig];
  TGraphAsymmErrors *a3CaloTrk_pt[n3Trig], *a3CaloTrk_eta[n3Trig];
  TGraphAsymmErrors *a4CaloTrk_pt[n4Trig], *a4CaloTrk_eta[n4Trig];
  TGraphAsymmErrors *a5CaloTrk_pt[n5Trig], *a5CaloTrk_eta[n5Trig];

  TGraphAsymmErrors *a2CaloTrk_MID_pt[n2Trig], *a2CaloTrk_MID_eta[n2Trig];
  TGraphAsymmErrors *a4CaloTrk_MID_pt[n4Trig], *a4CaloTrk_MID_eta[n4Trig];

  TGraphAsymmErrors *a2CaloGen_pt[n2Trig], *a2CaloGen_eta[n2Trig];
  TGraphAsymmErrors *a3CaloGen_pt[n3Trig], *a3CaloGen_eta[n3Trig];
  TGraphAsymmErrors *a4CaloGen_pt[n4Trig], *a4CaloGen_eta[n4Trig];
  TGraphAsymmErrors *a5CaloGen_pt[n5Trig], *a5CaloGen_eta[n5Trig];

  TGraphAsymmErrors *aPhoton_pt[nPhotonTrig], *aPhoton_eta[nPhotonTrig];
  TGraphAsymmErrors *aPhoton_MID_pt[nPhotonTrig], *aPhoton_MID_eta[nPhotonTrig];
  TGraphAsymmErrors *aDoublePhoton_pt[nPhotonTrig+2], *aDoublePhoton_eta[nPhotonTrig+2];

  for(int i = 0; i < n2Trig; ++i){
    a2Calo_pt[i] = new TGraphAsymmErrors();
    a2Calo_pt[i]->BayesDivide(hists2Calo_pt[i+1],hists2Calo_pt[0]);
    a2Calo_pt[i]->SetName(ak2CaloName[i]+"_asymm");
    a2Calo_eta[i] = new TGraphAsymmErrors();
    a2Calo_eta[i]->BayesDivide(hists2Calo_eta[i+1],hists2Calo_eta[0]);
    a2Calo_eta[i]->SetName(ak2CaloName[i]+"_eta_asymm");

    a2Calo_MID_pt[i] = new TGraphAsymmErrors();
    a2Calo_MID_pt[i]->BayesDivide(hists2Calo_MID_pt[i+1],hists2Calo_MID_pt[0]);
    a2Calo_MID_pt[i]->SetName(ak2CaloName_MID[i]+"_asymm");
    a2Calo_MID_eta[i] = new TGraphAsymmErrors();
    a2Calo_MID_eta[i]->BayesDivide(hists2Calo_MID_eta[i+1],hists2Calo_MID_eta[0]);
    a2Calo_MID_eta[i]->SetName(ak2CaloName_MID[i]+"_eta_asymm");

    a2CaloTrk_pt[i] = new TGraphAsymmErrors();
    a2CaloTrk_pt[i]->BayesDivide(hists2CaloTrk_pt[i+1],hists2CaloTrk_pt[0]);
    a2CaloTrk_pt[i]->SetName(ak2CaloTrkName[i]+"_asymm");
    a2CaloTrk_eta[i] = new TGraphAsymmErrors();
    a2CaloTrk_eta[i]->BayesDivide(hists2CaloTrk_eta[i+1],hists2CaloTrk_eta[0]);
    a2CaloTrk_eta[i]->SetName(ak2CaloTrkName[i]+"_eta_asymm");

    a2CaloTrk_MID_pt[i] = new TGraphAsymmErrors();
    a2CaloTrk_MID_pt[i]->BayesDivide(hists2CaloTrk_MID_pt[i+1],hists2CaloTrk_MID_pt[0]);
    a2CaloTrk_MID_pt[i]->SetName(ak2CaloTrkName_MID[i]+"_asymm");
    a2CaloTrk_MID_eta[i] = new TGraphAsymmErrors();
    a2CaloTrk_MID_eta[i]->BayesDivide(hists2CaloTrk_MID_eta[i+1],hists2CaloTrk_MID_eta[0]);
    a2CaloTrk_MID_eta[i]->SetName(ak2CaloTrkName_MID[i]+"_eta_asymm");

    a2CaloGen_pt[i] = new TGraphAsymmErrors();
    a2CaloGen_pt[i]->BayesDivide(hists2CaloGen_pt[i+1],hists2CaloGen_pt[0]);
    a2CaloGen_pt[i]->SetName(ak2CaloGenName[i]+"_asymm");
    a2CaloGen_eta[i] = new TGraphAsymmErrors();
    a2CaloGen_eta[i]->BayesDivide(hists2CaloGen_eta[i+1],hists2CaloGen_eta[0]);
    a2CaloGen_eta[i]->SetName(ak2CaloGenName[i]+"_eta_asymm");
  }

  for(int i = 0; i < n3Trig; ++i){
    a3Calo_pt[i] = new TGraphAsymmErrors();
    a3Calo_pt[i]->BayesDivide(hists3Calo_pt[i+1],hists3Calo_pt[0]);
    a3Calo_pt[i]->SetName(ak3CaloName[i]+"_asymm");
    a3Calo_eta[i] = new TGraphAsymmErrors();
    a3Calo_eta[i]->BayesDivide(hists3Calo_eta[i+1],hists3Calo_eta[0]);
    a3Calo_eta[i]->SetName(ak3CaloName[i]+"_eta_asymm");

    a3CaloTrk_pt[i] = new TGraphAsymmErrors();
    a3CaloTrk_pt[i]->BayesDivide(hists3CaloTrk_pt[i+1],hists3CaloTrk_pt[0]);
    a3CaloTrk_pt[i]->SetName(ak3CaloTrkName[i]+"_asymm");
    a3CaloTrk_eta[i] = new TGraphAsymmErrors();
    a3CaloTrk_eta[i]->BayesDivide(hists3CaloTrk_eta[i+1],hists3CaloTrk_eta[0]);
    a3CaloTrk_eta[i]->SetName(ak3CaloTrkName[i]+"_eta_asymm");

    a3CaloGen_pt[i] = new TGraphAsymmErrors();
    a3CaloGen_pt[i]->BayesDivide(hists3CaloGen_pt[i+1],hists3CaloGen_pt[0]);
    a3CaloGen_pt[i]->SetName(ak3CaloGenName[i]+"_asymm");
    a3CaloGen_eta[i] = new TGraphAsymmErrors();
    a3CaloGen_eta[i]->BayesDivide(hists3CaloGen_eta[i+1],hists3CaloGen_eta[0]);
    a3CaloGen_eta[i]->SetName(ak3CaloGenName[i]+"_eta_asymm");
  }

  for(int i = 0; i < n4Trig; ++i){
    a4Calo_pt[i] = new TGraphAsymmErrors();
    a4Calo_pt[i]->BayesDivide(hists4Calo_pt[i+1],hists4Calo_pt[0]);
    a4Calo_pt[i]->SetName(ak4CaloName[i]+"_asymm");
    a4Calo_eta[i] = new TGraphAsymmErrors();
    a4Calo_eta[i]->BayesDivide(hists4Calo_eta[i+1],hists4Calo_eta[0]);
    a4Calo_eta[i]->SetName(ak4CaloName[i]+"_eta_asymm");

    a4Calo_MID_pt[i] = new TGraphAsymmErrors();
    a4Calo_MID_pt[i]->BayesDivide(hists4Calo_MID_pt[i+1],hists4Calo_MID_pt[0]);
    a4Calo_MID_pt[i]->SetName(ak4CaloName_MID[i]+"_asymm");
    a4Calo_MID_eta[i] = new TGraphAsymmErrors();
    a4Calo_MID_eta[i]->BayesDivide(hists4Calo_MID_eta[i+1],hists4Calo_MID_eta[0]);
    a4Calo_MID_eta[i]->SetName(ak4CaloName_MID[i]+"_eta_asymm");

    a4CaloTrk_pt[i] = new TGraphAsymmErrors();
    a4CaloTrk_pt[i]->BayesDivide(hists4CaloTrk_pt[i+1],hists4CaloTrk_pt[0]);
    a4CaloTrk_pt[i]->SetName(ak4CaloTrkName[i]+"_asymm");
    a4CaloTrk_eta[i] = new TGraphAsymmErrors();
    a4CaloTrk_eta[i]->BayesDivide(hists4CaloTrk_eta[i+1],hists4CaloTrk_eta[0]);
    a4CaloTrk_eta[i]->SetName(ak4CaloTrkName[i]+"_eta_asymm");

    a4CaloTrk_MID_pt[i] = new TGraphAsymmErrors();
    a4CaloTrk_MID_pt[i]->BayesDivide(hists4CaloTrk_MID_pt[i+1],hists4CaloTrk_MID_pt[0]);
    a4CaloTrk_MID_pt[i]->SetName(ak4CaloTrkName_MID[i]+"_asymm");
    a4CaloTrk_MID_eta[i] = new TGraphAsymmErrors();
    a4CaloTrk_MID_eta[i]->BayesDivide(hists4CaloTrk_MID_eta[i+1],hists4CaloTrk_MID_eta[0]);
    a4CaloTrk_MID_eta[i]->SetName(ak4CaloTrkName_MID[i]+"_eta_asymm");

    a4CaloGen_pt[i] = new TGraphAsymmErrors();
    a4CaloGen_pt[i]->BayesDivide(hists4CaloGen_pt[i+1],hists4CaloGen_pt[0]);
    a4CaloGen_pt[i]->SetName(ak4CaloGenName[i]+"_asymm");
    a4CaloGen_eta[i] = new TGraphAsymmErrors();
    a4CaloGen_eta[i]->BayesDivide(hists4CaloGen_eta[i+1],hists4CaloGen_eta[0]);
    a4CaloGen_eta[i]->SetName(ak4CaloGenName[i]+"_eta_asymm");
  }

  a4Calo_3Jet_pt[0] = new TGraphAsymmErrors();
  a4Calo_3Jet_pt[0]->BayesDivide(hists4Calo_3Jet_pt[1],hists4Calo_3Jet_pt[0]);
  a4Calo_3Jet_pt[0]->SetName(ak4CaloName_3Jet[0]+"_asymm");
  a4Calo_3Jet_eta[0] = new TGraphAsymmErrors();
  a4Calo_3Jet_eta[0]->BayesDivide(hists4Calo_3Jet_eta[1],hists4Calo_3Jet_eta[0]);
  a4Calo_3Jet_eta[0]->SetName(ak4CaloName_3Jet[0]+"_eta_asymm");
  
  for(int i = 0; i < n5Trig; ++i){
    a5Calo_pt[i] = new TGraphAsymmErrors();
    a5Calo_pt[i]->BayesDivide(hists5Calo_pt[i+1],hists5Calo_pt[0]);
    a5Calo_pt[i]->SetName(ak5CaloName[i]+"_asymm");
    a5Calo_eta[i] = new TGraphAsymmErrors();
    a5Calo_eta[i]->BayesDivide(hists5Calo_eta[i+1],hists5Calo_eta[0]);
    a5Calo_eta[i]->SetName(ak5CaloName[i]+"_eta_asymm");

    a5CaloTrk_pt[i] = new TGraphAsymmErrors();
    a5CaloTrk_pt[i]->BayesDivide(hists5CaloTrk_pt[i+1],hists5CaloTrk_pt[0]);
    a5CaloTrk_pt[i]->SetName(ak5CaloTrkName[i]+"_asymm");
    a5CaloTrk_eta[i] = new TGraphAsymmErrors();
    a5CaloTrk_eta[i]->BayesDivide(hists5CaloTrk_eta[i+1],hists5CaloTrk_eta[0]);
    a5CaloTrk_eta[i]->SetName(ak5CaloTrkName[i]+"_eta_asymm");

    a5CaloGen_pt[i] = new TGraphAsymmErrors();
    a5CaloGen_pt[i]->BayesDivide(hists5CaloGen_pt[i+1],hists5CaloGen_pt[0]);
    a5CaloGen_pt[i]->SetName(ak5CaloGenName[i]+"_asymm");
    a5CaloGen_eta[i] = new TGraphAsymmErrors();
    a5CaloGen_eta[i]->BayesDivide(hists5CaloGen_eta[i+1],hists5CaloGen_eta[0]);
    a5CaloGen_eta[i]->SetName(ak5CaloGenName[i]+"_eta_asymm");
  }

  for(int i = 0; i < nPhotonTrig; ++i){
    aPhoton_pt[i] = new TGraphAsymmErrors();
    aPhoton_pt[i]->BayesDivide(histsPhoton_pt[i+1],histsPhoton_pt[0]);
    aPhoton_pt[i]->SetName(photonName[i]+"_asymm");
    aPhoton_eta[i] = new TGraphAsymmErrors();
    aPhoton_eta[i]->BayesDivide(histsPhoton_eta[i+1],histsPhoton_eta[0]);
    aPhoton_eta[i]->SetName(photonName[i]+"_eta_asymm");

    aPhoton_MID_pt[i] = new TGraphAsymmErrors();
    aPhoton_MID_pt[i]->BayesDivide(histsPhoton_MID_pt[i+1],histsPhoton_MID_pt[0]);
    aPhoton_MID_pt[i]->SetName(photonName_MID[i]+"_asymm");
    aPhoton_MID_eta[i] = new TGraphAsymmErrors();
    aPhoton_MID_eta[i]->BayesDivide(histsPhoton_MID_eta[i+1],histsPhoton_MID_eta[0]);
    aPhoton_MID_eta[i]->SetName(photonName_MID[i]+"_eta_asymm");
  }

  for(int i = 0; i < nPhotonTrig+2; ++i){
    aDoublePhoton_pt[i] = new TGraphAsymmErrors();
    aDoublePhoton_pt[i]->BayesDivide(histsDoublePhoton_pt[i+1],histsDoublePhoton_pt[0]);
    aDoublePhoton_pt[i]->SetName(photonDoubleName[i]+"_asymm");
    aDoublePhoton_eta[i] = new TGraphAsymmErrors();
    aDoublePhoton_eta[i]->BayesDivide(histsDoublePhoton_eta[i+1],histsDoublePhoton_eta[0]);
    aDoublePhoton_eta[i]->SetName(photonDoubleName[i]+"_eta_asymm");
  }

  //save output
  TFile *outFile = TFile::Open("jetTurnOn_HI.root","RECREATE");
  outFile->cd();

  hiBin_h->Write("", TObject::kOverwrite);
  for(Int_t iter = 0; iter < 4; iter++){
    rateHist_p[iter]->Write("", TObject::kOverwrite);
  }

  hists2Calo_pt[0]->Write("", TObject::kOverwrite);
  hists2Calo_eta[0]->Write("", TObject::kOverwrite);

  hists2Calo_MID_pt[0]->Write("", TObject::kOverwrite);
  hists2Calo_MID_eta[0]->Write("", TObject::kOverwrite);

  hists2CaloTrk_pt[0]->Write("", TObject::kOverwrite);
  hists2CaloTrk_eta[0]->Write("", TObject::kOverwrite);

  hists2CaloTrk_MID_pt[0]->Write("", TObject::kOverwrite);
  hists2CaloTrk_MID_eta[0]->Write("", TObject::kOverwrite);

  hists2CaloGen_pt[0]->Write("", TObject::kOverwrite);
  hists2CaloGen_eta[0]->Write("", TObject::kOverwrite);

  for(int i = 0; i < n2Trig; ++i){
    hists2Calo_pt[i+1]->Write("", TObject::kOverwrite);
    hists2Calo_eta[i+1]->Write("", TObject::kOverwrite);
    a2Calo_pt[i]->Write("", TObject::kOverwrite);
    a2Calo_eta[i]->Write("", TObject::kOverwrite);

    hists2Calo_MID_pt[i+1]->Write("", TObject::kOverwrite);
    hists2Calo_MID_eta[i+1]->Write("", TObject::kOverwrite);
    a2Calo_MID_pt[i]->Write("", TObject::kOverwrite);
    a2Calo_MID_eta[i]->Write("", TObject::kOverwrite);

    hists2CaloTrk_pt[i+1]->Write("", TObject::kOverwrite);
    hists2CaloTrk_eta[i+1]->Write("", TObject::kOverwrite);
    a2CaloTrk_pt[i]->Write("", TObject::kOverwrite);
    a2CaloTrk_eta[i]->Write("", TObject::kOverwrite);

    hists2CaloTrk_MID_pt[i+1]->Write("", TObject::kOverwrite);
    hists2CaloTrk_MID_eta[i+1]->Write("", TObject::kOverwrite);
    a2CaloTrk_MID_pt[i]->Write("", TObject::kOverwrite);
    a2CaloTrk_MID_eta[i]->Write("", TObject::kOverwrite);

    hists2CaloGen_pt[i+1]->Write("", TObject::kOverwrite);
    hists2CaloGen_eta[i+1]->Write("", TObject::kOverwrite);
    a2CaloGen_pt[i]->Write("", TObject::kOverwrite);
    a2CaloGen_eta[i]->Write("", TObject::kOverwrite);
  }

  hists3Calo_pt[0]->Write("", TObject::kOverwrite);
  hists3Calo_eta[0]->Write("", TObject::kOverwrite);

  hists3CaloTrk_pt[0]->Write("", TObject::kOverwrite);
  hists3CaloTrk_eta[0]->Write("", TObject::kOverwrite);

  hists3CaloGen_pt[0]->Write("", TObject::kOverwrite);
  hists3CaloGen_eta[0]->Write("", TObject::kOverwrite);

  for(int i = 0; i < n3Trig; ++i){
    hists3Calo_pt[i+1]->Write("", TObject::kOverwrite);
    hists3Calo_eta[i+1]->Write("", TObject::kOverwrite);
    a3Calo_pt[i]->Write("", TObject::kOverwrite);
    a3Calo_eta[i]->Write("", TObject::kOverwrite);

    hists3CaloTrk_pt[i+1]->Write("", TObject::kOverwrite);
    hists3CaloTrk_eta[i+1]->Write("", TObject::kOverwrite);
    a3CaloTrk_pt[i]->Write("", TObject::kOverwrite);
    a3CaloTrk_eta[i]->Write("", TObject::kOverwrite);

    hists3CaloGen_pt[i+1]->Write("", TObject::kOverwrite);
    hists3CaloGen_eta[i+1]->Write("", TObject::kOverwrite);
    a3CaloGen_pt[i]->Write("", TObject::kOverwrite);
    a3CaloGen_eta[i]->Write("", TObject::kOverwrite);
  }

  hists4Calo_pt[0]->Write("", TObject::kOverwrite);
  hists4Calo_eta[0]->Write("", TObject::kOverwrite);

  hists4Calo_MID_pt[0]->Write("", TObject::kOverwrite);
  hists4Calo_MID_eta[0]->Write("", TObject::kOverwrite);

  hists4Calo_3Jet_pt[0]->Write("", TObject::kOverwrite);
  hists4Calo_3Jet_eta[0]->Write("", TObject::kOverwrite);

  hists4CaloTrk_pt[0]->Write("", TObject::kOverwrite);
  hists4CaloTrk_eta[0]->Write("", TObject::kOverwrite);

  hists4CaloTrk_MID_pt[0]->Write("", TObject::kOverwrite);
  hists4CaloTrk_MID_eta[0]->Write("", TObject::kOverwrite);

  hists4CaloGen_pt[0]->Write("", TObject::kOverwrite);
  hists4CaloGen_eta[0]->Write("", TObject::kOverwrite);

  for(int i = 0; i < n4Trig; ++i){
    hists4Calo_pt[i+1]->Write("", TObject::kOverwrite);
    hists4Calo_eta[i+1]->Write("", TObject::kOverwrite);
    a4Calo_pt[i]->Write("", TObject::kOverwrite);
    a4Calo_eta[i]->Write("", TObject::kOverwrite);

    hists4Calo_MID_pt[i+1]->Write("", TObject::kOverwrite);
    hists4Calo_MID_eta[i+1]->Write("", TObject::kOverwrite);
    a4Calo_MID_pt[i]->Write("", TObject::kOverwrite);
    a4Calo_MID_eta[i]->Write("", TObject::kOverwrite);

    hists4CaloTrk_pt[i+1]->Write("", TObject::kOverwrite);
    hists4CaloTrk_eta[i+1]->Write("", TObject::kOverwrite);
    a4CaloTrk_pt[i]->Write("", TObject::kOverwrite);
    a4CaloTrk_eta[i]->Write("", TObject::kOverwrite);

    hists4CaloTrk_MID_pt[i+1]->Write("", TObject::kOverwrite);
    hists4CaloTrk_MID_eta[i+1]->Write("", TObject::kOverwrite);
    a4CaloTrk_MID_pt[i]->Write("", TObject::kOverwrite);
    a4CaloTrk_MID_eta[i]->Write("", TObject::kOverwrite);

    hists4CaloGen_pt[i+1]->Write("", TObject::kOverwrite);
    hists4CaloGen_eta[i+1]->Write("", TObject::kOverwrite);
    a4CaloGen_pt[i]->Write("", TObject::kOverwrite);
    a4CaloGen_eta[i]->Write("", TObject::kOverwrite);
  }

  hists4Calo_3Jet_pt[1]->Write("", TObject::kOverwrite);
  hists4Calo_3Jet_eta[1]->Write("", TObject::kOverwrite);
  a4Calo_3Jet_pt[0]->Write("", TObject::kOverwrite);
  a4Calo_3Jet_eta[0]->Write("", TObject::kOverwrite);
  

  hists5Calo_pt[0]->Write("", TObject::kOverwrite);
  hists5Calo_eta[0]->Write("", TObject::kOverwrite);

  hists5CaloTrk_pt[0]->Write("", TObject::kOverwrite);
  hists5CaloTrk_eta[0]->Write("", TObject::kOverwrite);

  hists5CaloGen_pt[0]->Write("", TObject::kOverwrite);
  hists5CaloGen_eta[0]->Write("", TObject::kOverwrite);

  for(int i = 0; i < n5Trig; ++i){
    hists5Calo_pt[i+1]->Write("", TObject::kOverwrite);
    hists5Calo_eta[i+1]->Write("", TObject::kOverwrite);
    a5Calo_pt[i]->Write("", TObject::kOverwrite);
    a5Calo_eta[i]->Write("", TObject::kOverwrite);

    hists5CaloTrk_pt[i+1]->Write("", TObject::kOverwrite);
    hists5CaloTrk_eta[i+1]->Write("", TObject::kOverwrite);
    a5CaloTrk_pt[i]->Write("", TObject::kOverwrite);
    a5CaloTrk_eta[i]->Write("", TObject::kOverwrite);

    hists5CaloGen_pt[i+1]->Write("", TObject::kOverwrite);
    hists5CaloGen_eta[i+1]->Write("", TObject::kOverwrite);
    a5CaloGen_pt[i]->Write("", TObject::kOverwrite);
    a5CaloGen_eta[i]->Write("", TObject::kOverwrite);
  }

  histsPhoton_pt[0]->Write("", TObject::kOverwrite);
  histsPhoton_eta[0]->Write("", TObject::kOverwrite);

  histsPhoton_MID_pt[0]->Write("", TObject::kOverwrite);
  histsPhoton_MID_eta[0]->Write("", TObject::kOverwrite);

  histsDoublePhoton_pt[0]->Write("", TObject::kOverwrite);
  histsDoublePhoton_eta[0]->Write("", TObject::kOverwrite);

  for(int i = 0; i < nPhotonTrig; ++i){
    histsPhoton_pt[i+1]->Write("", TObject::kOverwrite);
    histsPhoton_eta[i+1]->Write("", TObject::kOverwrite);
    aPhoton_pt[i]->Write("", TObject::kOverwrite);
    aPhoton_eta[i]->Write("", TObject::kOverwrite);

    histsPhoton_MID_pt[i+1]->Write("", TObject::kOverwrite);
    histsPhoton_MID_eta[i+1]->Write("", TObject::kOverwrite);
    aPhoton_MID_pt[i]->Write("", TObject::kOverwrite);
    aPhoton_MID_eta[i]->Write("", TObject::kOverwrite);
  }

  for(int i = 0; i < nPhotonTrig+2; ++i){
    histsDoublePhoton_pt[i+1]->Write("", TObject::kOverwrite);
    histsDoublePhoton_eta[i+1]->Write("", TObject::kOverwrite);
    aDoublePhoton_pt[i]->Write("", TObject::kOverwrite);
    aDoublePhoton_eta[i]->Write("", TObject::kOverwrite);
  }


  for(Int_t iter = 0; iter < 3; iter++){
    recoLeadJet3_p[iter]->Scale(1./matched);
    recoLeadJet3_p[iter]->Write("", TObject::kOverwrite);

    genLeadJet3_p[iter]->Scale(1./matched);
    genLeadJet3_p[iter]->Write("", TObject::kOverwrite);

    recoLeadJet4_p[iter]->Scale(1./matched);
    recoLeadJet4_p[iter]->Write("", TObject::kOverwrite);

    genLeadJet4_p[iter]->Scale(1./matched);
    genLeadJet4_p[iter]->Write("", TObject::kOverwrite);
  }

  TFile* out2 = new TFile("rateGen4Jet_TEMP.root", "UPDATE");

  for(Int_t iter = 0; iter < genLeadJet4_p[0]->GetNbinsX(); iter++){
    genLeadJet3_p[0]->SetBinContent(iter+1, genLeadJet3_p[0]->GetBinContent(iter+1)*30000.);
    recoLeadJet3_p[0]->SetBinContent(iter+1, recoLeadJet3_p[0]->GetBinContent(iter+1)*30000.);

    genLeadJet4_p[0]->SetBinContent(iter+1, genLeadJet4_p[0]->GetBinContent(iter+1)*30000.);
    recoLeadJet4_p[0]->SetBinContent(iter+1, recoLeadJet4_p[0]->GetBinContent(iter+1)*30000.);
  }

  genLeadJet3_p[0]->SetXTitle("Lead gen3Pt");
  genLeadJet3_p[0]->SetYTitle("Differential Rate (Hz)");
  genLeadJet3_p[0]->Write("gen3DiffRate_h", TObject::kOverwrite);

  recoLeadJet3_p[0]->SetXTitle("Lead reco3Pt");
  recoLeadJet3_p[0]->SetYTitle("Differential Rate (Hz)");
  recoLeadJet3_p[0]->Write("reco3DiffRate_h", TObject::kOverwrite);

  genLeadJet4_p[0]->SetXTitle("Lead gen4Pt");
  genLeadJet4_p[0]->SetYTitle("Differential Rate (Hz)");
  genLeadJet4_p[0]->Write("gen4DiffRate_h", TObject::kOverwrite);

  recoLeadJet4_p[0]->SetXTitle("Lead reco4Pt");
  recoLeadJet4_p[0]->SetYTitle("Differential Rate (Hz)");
  recoLeadJet4_p[0]->Write("reco4DiffRate_h", TObject::kOverwrite);

  for(Int_t iter = 0; iter < genLeadJet4_p[0]->GetNbinsX(); iter++){
    genLeadJet4_p[0]->SetBinContent(iter+1, genLeadJet4_p[0]->Integral(iter+2, 200));
    recoLeadJet4_p[0]->SetBinContent(iter+1, recoLeadJet4_p[0]->Integral(iter+2, 200));

    genLeadJet3_p[0]->SetBinContent(iter+1, genLeadJet3_p[0]->Integral(iter+2, 200));
    recoLeadJet3_p[0]->SetBinContent(iter+1, recoLeadJet3_p[0]->Integral(iter+2, 200));
  }

  genLeadJet3_p[0]->SetYTitle("Integrated Rate (Hz)");
  genLeadJet3_p[0]->Write("gen3IntegRate_h", TObject::kOverwrite);

  recoLeadJet3_p[0]->SetYTitle("Integrated Rate (Hz)");
  recoLeadJet3_p[0]->Write("reco3IntegRate_h", TObject::kOverwrite);

  genLeadJet4_p[0]->SetYTitle("Integrated Rate (Hz)");
  genLeadJet4_p[0]->Write("gen4IntegRate_h", TObject::kOverwrite);

  recoLeadJet4_p[0]->SetYTitle("Integrated Rate (Hz)");
  recoLeadJet4_p[0]->Write("reco4IntegRate_h", TObject::kOverwrite);

  for(Int_t iter = 0; iter < totIntPt80_h->GetNbinsX(); iter++){
    totIntPt80_h->SetBinContent(iter+1, totIntPt80_h->Integral(iter+2, totIntPt80_h->GetNbinsX()));
    totIntPt80_h->SetBinError(iter+1, TMath::Sqrt(totIntPt80_h->Integral(iter+2, totIntPt80_h->GetNbinsX())));

    maxMissedIntPt80_h->SetBinContent(iter+1, maxMissedIntPt80_h->Integral(iter+2, maxMissedIntPt80_h->GetNbinsX()));
    maxMissedIntPt80_h->SetBinError(iter+1, TMath::Sqrt(maxMissedIntPt80_h->Integral(iter+2, maxMissedIntPt80_h->GetNbinsX())));
  }

  for(Int_t iter = 0; iter < totIntPt100_h->GetNbinsX(); iter++){
    totIntPt100_h->SetBinContent(iter+1, totIntPt100_h->Integral(iter+2, totIntPt100_h->GetNbinsX()));
    totIntPt100_h->SetBinError(iter+1, TMath::Sqrt(totIntPt100_h->Integral(iter+2, totIntPt100_h->GetNbinsX())));

    maxMissedIntPt100_h->SetBinContent(iter+1, maxMissedIntPt100_h->Integral(iter+2, maxMissedIntPt100_h->GetNbinsX()));
    maxMissedIntPt100_h->SetBinError(iter+1, TMath::Sqrt(maxMissedIntPt100_h->Integral(iter+2, maxMissedIntPt100_h->GetNbinsX())));
  }

  for(Int_t iter = 0; iter < totIntPt120_h->GetNbinsX(); iter++){
    totIntPt120_h->SetBinContent(iter+1, totIntPt120_h->Integral(iter+2, totIntPt120_h->GetNbinsX()));
    totIntPt120_h->SetBinError(iter+1, TMath::Sqrt(totIntPt120_h->Integral(iter+2, totIntPt120_h->GetNbinsX())));

    maxMissedIntPt120_h->SetBinContent(iter+1, maxMissedIntPt120_h->Integral(iter+2, maxMissedIntPt120_h->GetNbinsX()));
    maxMissedIntPt120_h->SetBinError(iter+1, TMath::Sqrt(maxMissedIntPt120_h->Integral(iter+2, maxMissedIntPt120_h->GetNbinsX())));
  }


  maxMissedPt80_h->Divide(totPt80_h);
  maxMissedEta80_h->Divide(totEta80_h);

  maxMissedPt100_h->Divide(totPt100_h);
  maxMissedEta100_h->Divide(totEta100_h);

  maxMissedPt120_h->Divide(totPt120_h);
  maxMissedEta120_h->Divide(totEta120_h);

  maxMissedIntPt80_h->Divide(totIntPt80_h);
  maxMissedIntPt100_h->Divide(totIntPt100_h);
  maxMissedIntPt120_h->Divide(totIntPt120_h);

  maxMissedHiBin80_h->Divide(totHiBin80_h);
  maxMissedHiBin100_h->Divide(totHiBin100_h);
  maxMissedHiBin120_h->Divide(totHiBin120_h);


  maxMissedPt80_h->Write("", TObject::kOverwrite);
  maxMissedEta80_h->Write("", TObject::kOverwrite);

  maxMissedPt100_h->GetXaxis()->CenterTitle();
  maxMissedPt100_h->SetXTitle("Offline Jet p_{T}");

  maxMissedPt100_h->GetYaxis()->CenterTitle();
  maxMissedPt100_h->SetYTitle("Missed Fraction (Diff.)");

  maxMissedPt100_h->Write("", TObject::kOverwrite);
  maxMissedEta100_h->Write("", TObject::kOverwrite);

  maxMissedPt120_h->Write("", TObject::kOverwrite);
  maxMissedEta120_h->Write("", TObject::kOverwrite);

  maxMissedIntPt100_h->GetXaxis()->CenterTitle();
  maxMissedIntPt100_h->SetXTitle("Offline Jet p_{T}");

  maxMissedIntPt100_h->GetYaxis()->CenterTitle();
  maxMissedIntPt100_h->SetYTitle("Missed Fraction (Int.)");

  maxMissedIntPt80_h->Write("", TObject::kOverwrite);
  maxMissedIntPt100_h->Write("", TObject::kOverwrite);
  maxMissedIntPt120_h->Write("", TObject::kOverwrite);

  maxMissedHiBin80_h->Write("", TObject::kOverwrite);
  maxMissedHiBin100_h->Write("", TObject::kOverwrite);
  maxMissedHiBin120_h->Write("", TObject::kOverwrite);


  out2->Close();
  delete out2;
  delete hiBin_h;

  HLTFile->Close();
  AnaFile->Close();
  outFile->Close();
}

int main()
{
  matchJetTree();
  return 0;
}
