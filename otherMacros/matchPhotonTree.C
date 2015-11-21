#include "EventMatchingCMS.h"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <iostream>
#include <TStyle.h>

const TString AnaFilename = "SimplePhotonAnalyzer_NeutrinoGun_20150326.root";
const TString AnaPhotonTreename = "SimplePhotonAnalyzer/PhotonTree";
//const TString HLTFilename = "hltbits_20150304_740pre7_RelValProdTTBar.root";
const TString HLTFilename = "hltbits_LowPU_20150603.root";

const int nBins = 100;
const double maxpt = 100;

void matchJetTree()
{
  gStyle->SetOptStat(0);

  TFile *HLTFile = TFile::Open(HLTFilename);
  TTree *HLTTree = (TTree*)HLTFile->Get("HltTree");

  ULong64_t hlt_event;
  Int_t hlt_run, hlt_lumi;

  const Int_t nL1Trig = 9;
  const Int_t nHLTTrig = 5;

  Int_t L1_SingleEG2 = 0;
  Int_t L1_SingleEG5 = 0;
  Int_t L1_SingleEG10 = 0;
  Int_t L1_SingleEG15 = 0;
  Int_t L1_SingleEG20 = 0;
  Int_t L1_SingleEG25 = 0;
  Int_t L1_SingleEG30 = 0;
  Int_t L1_SingleEG35 = 0;
  Int_t L1_SingleEG40 = 0;

  Int_t nL12 = 0;
  Int_t nL15 = 0;
  Int_t nL110 = 0;
  Int_t nL115 = 0;
  Int_t nL120 = 0;
  Int_t nL125 = 0;
  Int_t nL130 = 0;
  Int_t nL135 = 0;
  Int_t nL140 = 0;

  Int_t HLT_HISinglePhoton10_v1;
  Int_t HLT_HISinglePhoton15_v1;
  Int_t HLT_HISinglePhoton20_v1;
  Int_t HLT_HISinglePhoton40_v1;
  Int_t HLT_HISinglePhoton60_v1;

  Int_t nPhoton10 = 0;
  Int_t nPhoton15 = 0;
  Int_t nPhoton20 = 0;
  Int_t nPhoton40 = 0;
  Int_t nPhoton60 = 0;

  Int_t nFakePhoton10 = 0;
  Int_t nFakePhoton15 = 0;
  Int_t nFakePhoton20 = 0;
  Int_t nFakePhoton40 = 0;
  Int_t nFakePhoton60 = 0;

  TString L1Name[nL1Trig] = {"L1_SingleEG2", "L1_SingleEG5", "L1_SingleEG10", "L1_SingleEG15", "L1_SingleEG20", "L1_SingleEG25", "L1_SingleEG30", "L1_SingleEG35", "L1_SingleEG40"};

  TString photonName[nHLTTrig] = {"HLT_HISinglePhoton10_v2", "HLT_HISinglePhoton15_v2", "HLT_HISinglePhoton20_v2", "HLT_HISinglePhoton40_v2", "HLT_HISinglePhoton60_v2"};

  HLTTree->SetBranchStatus("*", 0);

  HLTTree->SetBranchStatus("Event", 1);
  HLTTree->SetBranchStatus("Run", 1);
  HLTTree->SetBranchStatus("LumiBlock", 1);

  HLTTree->SetBranchStatus(L1Name[0], 1);
  HLTTree->SetBranchStatus(L1Name[1], 1);
  HLTTree->SetBranchStatus(L1Name[2], 1);
  HLTTree->SetBranchStatus(L1Name[3], 1);
  HLTTree->SetBranchStatus(L1Name[4], 1);
  HLTTree->SetBranchStatus(L1Name[5], 1);
  HLTTree->SetBranchStatus(L1Name[6], 1);
  HLTTree->SetBranchStatus(L1Name[7], 1);
  HLTTree->SetBranchStatus(L1Name[8], 1);

  HLTTree->SetBranchStatus(photonName[0], 1);
  HLTTree->SetBranchStatus(photonName[1], 1);
  HLTTree->SetBranchStatus(photonName[2], 1);
  HLTTree->SetBranchStatus(photonName[3], 1);
  HLTTree->SetBranchStatus(photonName[4], 1);



  HLTTree->SetBranchAddress("Event",&hlt_event);
  HLTTree->SetBranchAddress("Run",&hlt_run);
  HLTTree->SetBranchAddress("LumiBlock",&hlt_lumi);

  HLTTree->SetBranchAddress(L1Name[0], &L1_SingleEG2);
  HLTTree->SetBranchAddress(L1Name[1], &L1_SingleEG5);
  HLTTree->SetBranchAddress(L1Name[2], &L1_SingleEG10);
  HLTTree->SetBranchAddress(L1Name[3], &L1_SingleEG15);
  HLTTree->SetBranchAddress(L1Name[4], &L1_SingleEG20);
  HLTTree->SetBranchAddress(L1Name[5], &L1_SingleEG25);
  HLTTree->SetBranchAddress(L1Name[6], &L1_SingleEG30);
  HLTTree->SetBranchAddress(L1Name[7], &L1_SingleEG35);
  HLTTree->SetBranchAddress(L1Name[8], &L1_SingleEG40);

  HLTTree->SetBranchAddress(photonName[0], &HLT_HISinglePhoton10_v1);
  HLTTree->SetBranchAddress(photonName[1], &HLT_HISinglePhoton15_v1);
  HLTTree->SetBranchAddress(photonName[2], &HLT_HISinglePhoton20_v1);
  HLTTree->SetBranchAddress(photonName[3], &HLT_HISinglePhoton40_v1);
  HLTTree->SetBranchAddress(photonName[4], &HLT_HISinglePhoton60_v1);


  TFile *AnaFile = TFile::Open(AnaFilename);
  TTree *AnaPhotonTree = (TTree*)AnaFile->Get(AnaPhotonTreename); // other option is SimplePhotonAnalyzer

  Int_t ana_event, ana_lumi;//, ana_run, ana_lumi;
  Int_t nPhotons;
  Double_t pt[500], eta[500], phi[500];

  AnaPhotonTree->SetBranchAddress("event", &ana_event);
  //  AnaPhotonTree->SetBranchAddress("run", &ana_run);
  AnaPhotonTree->SetBranchAddress("lumi", &ana_lumi);

  AnaPhotonTree->SetBranchAddress("nPhotons", &nPhotons);
  AnaPhotonTree->SetBranchAddress("pt", pt);
  AnaPhotonTree->SetBranchAddress("eta", eta);
  AnaPhotonTree->SetBranchAddress("phi", phi);

  //book histos
  TH1D *histsL1_pt[nL1Trig + 1], *histsL1_eta[nL1Trig + 1];
  TH1D *histsHLT_pt[nHLTTrig + 1], *histsHLT_eta[nHLTTrig + 1];
 
  histsL1_pt[0] = new TH1D("leadingL1_pt",";p_{T}^{jt}",nBins,0,maxpt);
  histsL1_eta[0] = new TH1D("leadingL1_eta",";#eta^{jt}",nBins,-5,5);

  histsHLT_pt[0] = new TH1D("leadingHLT_pt",";p_{T}^{jt}",nBins,0,maxpt);
  histsHLT_eta[0] = new TH1D("leadingHLT_eta",";#eta^{jt}",nBins,-5,5);

  for(int i = 0; i < nL1Trig; ++i){
    histsL1_pt[i+1] = (TH1D*)histsL1_pt[0]->Clone(L1Name[i]);
    histsL1_eta[i+1] = (TH1D*)histsL1_eta[0]->Clone(L1Name[i]+"eta");
  }

  for(int i = 0; i < nHLTTrig; ++i){
    histsHLT_pt[i+1] = (TH1D*)histsHLT_pt[0]->Clone(photonName[i]);
    histsHLT_eta[i+1] = (TH1D*)histsHLT_eta[0]->Clone(photonName[i]+"eta");
  }

  std::cout << "Events in HLT file: " << HLTTree->GetEntries() << std::endl;
  std::cout << "Events in Ana file: " << AnaPhotonTree->GetEntries() << std::endl;

  //make map
  EventMatchingCMS *matcher = new EventMatchingCMS();
  for(Long64_t entry = 0; entry < HLTTree->GetEntries(); ++entry)
  {
    HLTTree->GetEntry(entry);
    matcher->addEvent(hlt_event, hlt_lumi, 0, entry);
  }

  // analysis loop
  int matched = 0;
  for(Long64_t entry = 0; entry < AnaPhotonTree->GetEntries(); ++entry)
  {
    if(entry%100000 == 0) std::cout << entry << std::endl;

    AnaPhotonTree->GetEntry(entry);

    long long hlt_entry = matcher->retrieveEvent(ana_event, ana_lumi, 0);
    if(hlt_entry == -1) continue;
    HLTTree->GetEntry(hlt_entry);
    matched++;

    if(L1_SingleEG2) nL12++;
    if(L1_SingleEG5) nL15++;
    if(L1_SingleEG10) nL110++;
    if(L1_SingleEG15) nL115++;
    if(L1_SingleEG20) nL120++;
    if(L1_SingleEG25) nL125++;
    if(L1_SingleEG30) nL130++;
    if(L1_SingleEG35) nL135++;
    if(L1_SingleEG40) nL140++;

    if(HLT_HISinglePhoton10_v1) nPhoton10++;
    if(HLT_HISinglePhoton15_v1) nPhoton15++;
    if(HLT_HISinglePhoton20_v1) nPhoton20++;
    if(HLT_HISinglePhoton40_v1) nPhoton40++;
    if(HLT_HISinglePhoton60_v1) nPhoton60++;

    Double_t maxPhotonPt = -1;
    Double_t maxPhotonEta = -100;
    for(int i = 0; i < nPhotons; ++i){
      if(fabs(eta[i]) > 2.0) continue;
      if(pt[i] > maxPhotonPt){
	maxPhotonPt = pt[i];
	maxPhotonEta = eta[i];
      }
    }


    if(maxPhotonPt != -1){
      histsL1_pt[0]->Fill(maxPhotonPt);
      histsL1_eta[0]->Fill(maxPhotonEta);

      if(L1_SingleEG2){
	histsL1_pt[1]->Fill(maxPhotonPt);
	histsL1_eta[1]->Fill(maxPhotonEta);
      }
      if(L1_SingleEG5){
	histsL1_pt[2]->Fill(maxPhotonPt);
	histsL1_eta[2]->Fill(maxPhotonEta);
      }
      if(L1_SingleEG10){
	histsL1_pt[3]->Fill(maxPhotonPt);
	histsL1_eta[3]->Fill(maxPhotonEta);
      }
      if(L1_SingleEG15){
	histsL1_pt[4]->Fill(maxPhotonPt);
	histsL1_eta[4]->Fill(maxPhotonEta);
      }
      if(L1_SingleEG20){
	histsL1_pt[5]->Fill(maxPhotonPt);
	histsL1_eta[5]->Fill(maxPhotonEta);
      }
      if(L1_SingleEG25){
	histsL1_pt[6]->Fill(maxPhotonPt);
	histsL1_eta[6]->Fill(maxPhotonEta);
      }
      if(L1_SingleEG30){
	histsL1_pt[7]->Fill(maxPhotonPt);
	histsL1_eta[7]->Fill(maxPhotonEta);
      }
      if(L1_SingleEG35){
	histsL1_pt[8]->Fill(maxPhotonPt);
	histsL1_eta[8]->Fill(maxPhotonEta);
      }
      if(L1_SingleEG40){
	histsL1_pt[9]->Fill(maxPhotonPt);
	histsL1_eta[9]->Fill(maxPhotonEta);
      }

      histsHLT_pt[0]->Fill(maxPhotonPt);
      histsHLT_eta[0]->Fill(maxPhotonEta);

      if(HLT_HISinglePhoton10_v1){
	histsHLT_pt[1]->Fill(maxPhotonPt);
	histsHLT_eta[1]->Fill(maxPhotonEta);
      }
      if(HLT_HISinglePhoton15_v1){
	histsHLT_pt[2]->Fill(maxPhotonPt);
	histsHLT_eta[2]->Fill(maxPhotonEta);
      }
      if(HLT_HISinglePhoton20_v1){
	histsHLT_pt[3]->Fill(maxPhotonPt);
	histsHLT_eta[3]->Fill(maxPhotonEta);
      }
      if(HLT_HISinglePhoton40_v1){
	histsHLT_pt[4]->Fill(maxPhotonPt);
	histsHLT_eta[4]->Fill(maxPhotonEta);
      }
      if(HLT_HISinglePhoton60_v1){
	histsHLT_pt[5]->Fill(maxPhotonPt);
	histsHLT_eta[5]->Fill(maxPhotonEta);
      }    
    }  
    else{
      if(HLT_HISinglePhoton10_v1) nFakePhoton10++;
      if(HLT_HISinglePhoton15_v1) nFakePhoton15++;
      if(HLT_HISinglePhoton20_v1) nFakePhoton20++;
      if(HLT_HISinglePhoton40_v1) nFakePhoton40++;
      if(HLT_HISinglePhoton60_v1) nFakePhoton60++;
    }
  }

  std::cout << "Events matched: " << matched << std::endl;

  std::cout << "True trigger rates: " << std::endl;

  std::cout << "L1_SingleEG2: " << nL12 << std::endl;
  std::cout << "L1_SingleEG5: " << nL15 << std::endl;
  std::cout << "L1_SingleEG10: " << nL110 << std::endl;
  std::cout << "L1_SingleEG15: " << nL115 << std::endl;
  std::cout << "L1_SingleEG20: " << nL120 << std::endl;
  std::cout << "L1_SingleEG25: " << nL125 << std::endl;
  std::cout << "L1_SingleEG30: " << nL130 << std::endl;
  std::cout << "L1_SingleEG35: " << nL135 << std::endl;
  std::cout << "L1_SingleEG40: " << nL140 << std::endl;

  std::cout << "HLT_HISinglePhoton10_v1: " << nPhoton10 << ", " << nFakePhoton10 << std::endl;
  std::cout << "HLT_HISinglePhoton15_v1: " << nPhoton15 << ", " << nFakePhoton15 << std::endl;
  std::cout << "HLT_HISinglePhoton20_v1: " << nPhoton20 << ", " << nFakePhoton20 << std::endl;
  std::cout << "HLT_HISinglePhoton40_v1: " << nPhoton40 << ", " << nFakePhoton40 << std::endl;
  std::cout << "HLT_HISinglePhoton60_v1: " << nPhoton60 << ", " << nFakePhoton60 << std::endl;

  std::cout << std::endl;
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
  TGraphAsymmErrors *aL1_pt[nL1Trig], *aL1_eta[nL1Trig];
  TGraphAsymmErrors *aHLT_pt[nHLTTrig], *aHLT_eta[nHLTTrig];

  for(int i = 0; i < nL1Trig; ++i){
    aL1_pt[i] = new TGraphAsymmErrors();
    aL1_pt[i]->BayesDivide(histsL1_pt[i+1],histsL1_pt[0]);
    aL1_pt[i]->SetName(L1Name[i]+"_asymm");
    aL1_eta[i] = new TGraphAsymmErrors();
    aL1_eta[i]->BayesDivide(histsL1_eta[i+1],histsL1_eta[0]);
    aL1_eta[i]->SetName(L1Name[i]+"_eta_asymm");
  }

  for(int i = 0; i < nHLTTrig; ++i){
    aHLT_pt[i] = new TGraphAsymmErrors();
    aHLT_pt[i]->BayesDivide(histsHLT_pt[i+1],histsHLT_pt[0]);
    aHLT_pt[i]->SetName(photonName[i]+"_asymm");
    aHLT_eta[i] = new TGraphAsymmErrors();
    aHLT_eta[i]->BayesDivide(histsHLT_eta[i+1],histsHLT_eta[0]);
    aHLT_eta[i]->SetName(photonName[i]+"_eta_asymm");
  }

  //save output
  TFile *outFile = TFile::Open("photonTurnOn.root","RECREATE");
  outFile->cd();

  histsL1_pt[0]->Write("", TObject::kOverwrite);
  histsL1_eta[0]->Write("", TObject::kOverwrite);

  histsHLT_pt[0]->Write("", TObject::kOverwrite);
  histsHLT_eta[0]->Write("", TObject::kOverwrite);

  for(int i = 0; i < nL1Trig; ++i){
    histsL1_pt[i+1]->Write("", TObject::kOverwrite);
    histsL1_eta[i+1]->Write("", TObject::kOverwrite);
    aL1_pt[i]->Write("", TObject::kOverwrite);
    aL1_eta[i]->Write("", TObject::kOverwrite);
  }

  for(int i = 0; i < nHLTTrig; ++i){
    histsHLT_pt[i+1]->Write("", TObject::kOverwrite);
    histsHLT_eta[i+1]->Write("", TObject::kOverwrite);
    aHLT_pt[i]->Write("", TObject::kOverwrite);
    aHLT_eta[i]->Write("", TObject::kOverwrite);
  }

  HLTFile->Close();
  AnaFile->Close();
  outFile->Close();
}

int main()
{
  matchJetTree();
  return 0;
}
