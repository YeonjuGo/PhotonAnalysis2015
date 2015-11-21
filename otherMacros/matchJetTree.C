#include "EventMatchingCMS.h"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <iostream>
#include <TStyle.h>
#include "TMath.h"

//const TString AnaFilename = "SimpleJetAnalyzer_20150304_740pre7_RelValProdTTBar.root";
const TString AnaFilename = "SimpleJetAnalyzer_NeutrinoGun_20150302.root";
const TString AnaCaloTreename = "ak5CaloJetAnalyzer/t";
const TString AnaPFTreename = "ak5PFJetAnalyzer/t";
//const TString HLTFilename = "hltbits_20150304_740pre7_RelValProdTTBar.root";
const TString HLTFilename = "hltbits_LowPU_20150616_SEED10_MBSEED.root";

const int nBins = 150;
const double maxpt = 150;

void matchJetTree()
{
  gStyle->SetOptStat(0);

  TFile *HLTFile = TFile::Open(HLTFilename);
  TTree *HLTTree = (TTree*)HLTFile->Get("HltTree");

  ULong64_t hlt_event;
  Int_t hlt_run, hlt_lumi;

  const Int_t nPtTrigCalo = 5;
  const Int_t nPtTrigPF = 4;

  Int_t HLT_AK4CaloJet30_v2;
  Int_t HLT_AK4CaloJet40_v2;
  Int_t HLT_AK4CaloJet50_v2;
  Int_t HLT_AK4CaloJet80_v2;
  Int_t HLT_AK4CaloJet100_v2;

  Int_t HLT_AK4PFJet30_v2;
  Int_t HLT_AK4PFJet50_v2;
  Int_t HLT_AK4PFJet80_v2;
  Int_t HLT_AK4PFJet100_v2;

  Int_t n4CaloJet30 = 0;
  Int_t n4CaloJet40 = 0;
  Int_t n4CaloJet50 = 0;
  Int_t n4CaloJet80 = 0;
  Int_t n4CaloJet100 = 0;

  Int_t n4PFJet30 = 0;
  Int_t n4PFJet50 = 0;
  Int_t n4PFJet80 = 0;
  Int_t n4PFJet100 = 0;

  Int_t nFake4CaloJet30 = 0;
  Int_t nFake4CaloJet40 = 0;
  Int_t nFake4CaloJet50 = 0;
  Int_t nFake4CaloJet80 = 0;
  Int_t nFake4CaloJet100 = 0;

  Int_t nFake4PFJet30 = 0;
  Int_t nFake4PFJet50 = 0;
  Int_t nFake4PFJet80 = 0;
  Int_t nFake4PFJet100 = 0;

  TString ak4CaloName[nPtTrigCalo] = {"HLT_AK4CaloJet30_v2", "HLT_AK4CaloJet40_v2", "HLT_AK4CaloJet50_v2", "HLT_AK4CaloJet80_v2", "HLT_AK4CaloJet100_v2"};

  TString ak4PFName[nPtTrigPF] = {"HLT_AK4PFJet30_v2", "HLT_AK4PFJet50_v2", "HLT_AK4PFJet80_v2", "HLT_AK4PFJet100_v2"};

  HLTTree->SetBranchAddress("Event",&hlt_event);
  HLTTree->SetBranchAddress("Run",&hlt_run);
  HLTTree->SetBranchAddress("LumiBlock",&hlt_lumi);

  HLTTree->SetBranchAddress(ak4CaloName[0], &HLT_AK4CaloJet30_v2);
  HLTTree->SetBranchAddress(ak4CaloName[1], &HLT_AK4CaloJet40_v2);
  HLTTree->SetBranchAddress(ak4CaloName[2], &HLT_AK4CaloJet50_v2);
  HLTTree->SetBranchAddress(ak4CaloName[3], &HLT_AK4CaloJet80_v2);
  HLTTree->SetBranchAddress(ak4CaloName[4], &HLT_AK4CaloJet100_v2);

  HLTTree->SetBranchAddress(ak4PFName[0], &HLT_AK4PFJet30_v2);
  HLTTree->SetBranchAddress(ak4PFName[1], &HLT_AK4PFJet50_v2);
  HLTTree->SetBranchAddress(ak4PFName[2], &HLT_AK4PFJet80_v2);
  HLTTree->SetBranchAddress(ak4PFName[3], &HLT_AK4PFJet100_v2);


  TFile *AnaFile = TFile::Open(AnaFilename);
  TTree *AnaCaloTree = (TTree*)AnaFile->Get(AnaCaloTreename); // other option is SimplePhotonAnalyzer
  TTree *AnaPFTree = (TTree*)AnaFile->Get(AnaPFTreename); // other option is SimplePhotonAnalyzer

  Int_t ana_event;//, ana_run, ana_lumi;
  Int_t nCaloref;
  Float_t jtCalopt[500], jtCaloeta[500], jtCalophi[500];

  Int_t ngen;
  Float_t genpt[500], geneta[500], genphi[500];

  AnaCaloTree->SetBranchAddress("evt", &ana_event);
  //  AnaCaloTree->SetBranchAddress("run", &ana_run);
  //  AnaCaloTree->SetBranchAddress("lumi", &ana_lumi);

  AnaCaloTree->SetBranchAddress("nref", &nCaloref);
  AnaCaloTree->SetBranchAddress("jtpt", jtCalopt);
  AnaCaloTree->SetBranchAddress("jteta", jtCaloeta);
  AnaCaloTree->SetBranchAddress("jtphi", jtCalophi);

  AnaCaloTree->SetBranchAddress("ngen", &ngen);
  AnaCaloTree->SetBranchAddress("genpt", genpt);
  AnaCaloTree->SetBranchAddress("geneta", geneta);
  AnaCaloTree->SetBranchAddress("genphi", genphi);

  Int_t nPFref;
  Int_t pf_event;
  Float_t jtPFpt[500], jtPFeta[500], jtPFphi[500];

  AnaPFTree->SetBranchAddress("evt", &pf_event);

  AnaPFTree->SetBranchAddress("nref", &nPFref);
  AnaPFTree->SetBranchAddress("jtpt", jtPFpt);
  AnaPFTree->SetBranchAddress("jteta", jtPFeta);
  AnaPFTree->SetBranchAddress("jtphi", jtPFphi);

  //book histos
  TH1D *hists4Calo_pt[nPtTrigCalo + 1], *hists4Calo_eta[nPtTrigCalo + 1];
  TH1D *hists4PF_pt[nPtTrigPF + 1], *hists4PF_eta[nPtTrigPF + 1];

  hists4Calo_pt[0] = new TH1D("leading_ak4Calo_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists4Calo_eta[0] = new TH1D("leading_ak4Calo_eta",";#eta^{jt}",nBins,-5,5);
  hists4PF_pt[0] = new TH1D("leading_ak4PF_pt",";p_{T}^{jt}",nBins,0,maxpt);
  hists4PF_eta[0] = new TH1D("leading_ak4PF_eta",";#eta^{jt}",nBins,-5,5);

  for(int i = 0; i < nPtTrigCalo; ++i){
    hists4Calo_pt[i+1] = (TH1D*)hists4Calo_pt[0]->Clone(ak4CaloName[i]);
    hists4Calo_eta[i+1] = (TH1D*)hists4Calo_eta[0]->Clone(ak4CaloName[i]+"eta");
  }

  for(int i = 0; i < nPtTrigPF; ++i){
    hists4PF_pt[i+1] = (TH1D*)hists4PF_pt[0]->Clone(ak4PFName[i]);
    hists4PF_eta[i+1] = (TH1D*)hists4PF_eta[0]->Clone(ak4PFName[i]+"eta");
  }

  std::cout << "Events in HLT file: " << HLTTree->GetEntries() << std::endl;
  std::cout << "Events in Ana file: " << AnaCaloTree->GetEntries() << std::endl;

  //make map
  EventMatchingCMS *matcher = new EventMatchingCMS();
  for(Long64_t entry = 0; entry < HLTTree->GetEntries(); ++entry)
  {
    HLTTree->GetEntry(entry);
    matcher->addEvent(hlt_event, 0, 0, entry);
  }

  // analysis loop
  int matched = 0;
  for(Long64_t entry = 0; entry < AnaCaloTree->GetEntries(); ++entry)
  {
    if(entry%100000 == 0) std::cout << entry << std::endl;

    AnaCaloTree->GetEntry(entry);
    AnaPFTree->GetEntry(entry);

    if(ana_event != pf_event) std::cout << "MOTHERFUCKER THERE IS A PROBLEM" << std::endl;

    long long hlt_entry = matcher->retrieveEvent(ana_event, 0, 0);
    if(hlt_entry == -1) continue;
    HLTTree->GetEntry(hlt_entry);
    matched++;

    if(HLT_AK4CaloJet30_v2) n4CaloJet30++;
    if(HLT_AK4CaloJet40_v2) n4CaloJet40++;
    if(HLT_AK4CaloJet50_v2) n4CaloJet50++;
    if(HLT_AK4CaloJet80_v2) n4CaloJet80++;
    if(HLT_AK4CaloJet100_v2) n4CaloJet100++;

    if(HLT_AK4PFJet30_v2) n4PFJet30++;
    if(HLT_AK4PFJet50_v2) n4PFJet50++;
    if(HLT_AK4PFJet80_v2) n4PFJet80++;
    if(HLT_AK4PFJet100_v2) n4PFJet100++;

    Double_t maxCaloAnaPt = -1;
    Double_t maxCaloAnaEta = -100;
    for(int i = 0; i < nCaloref; ++i)
    {
      if(fabs(jtCaloeta[i]) > 2.0) continue;
      if(jtCalopt[i] > maxCaloAnaPt)
      {
	maxCaloAnaPt = jtCalopt[i];
	maxCaloAnaEta = jtCaloeta[i];
      }
    }

    Double_t maxPFAnaPt = -1;
    Double_t maxPFAnaEta = -100;
    for(int i = 0; i < nPFref; ++i)
    {
      if(fabs(jtPFeta[i]) > 2.0) continue;
      if(jtPFpt[i] > maxPFAnaPt)
      {
	maxPFAnaPt = jtPFpt[i];
	maxPFAnaEta = jtPFeta[i];
      }
    }

    Double_t maxGenPt = -1;
    //    Double_t maxGenEta = -100;
    for(int i = 0; i < ngen; ++i){
      if(fabs(geneta[i]) > 2.0) continue;
      if(genpt[i] > maxGenPt && genpt[i] > 1){
	maxGenPt = genpt[i];
	//	maxGenEta = geneta[i];
      }
    }

    //    maxCaloAnaPt = maxGenPt;
    //    maxCaloAnaEta = maxGenEta;

    //    maxPFAnaPt = maxGenPt;
    //    maxPFAnaEta = maxGenEta;

    if(maxCaloAnaPt != -1){
      hists4Calo_pt[0]->Fill(maxCaloAnaPt);
      hists4Calo_eta[0]->Fill(maxCaloAnaEta);

      if(HLT_AK4CaloJet30_v2){
	hists4Calo_pt[1]->Fill(maxCaloAnaPt);
	hists4Calo_eta[1]->Fill(maxCaloAnaEta);
      }
      if(HLT_AK4CaloJet40_v2){
	hists4Calo_pt[2]->Fill(maxCaloAnaPt);
	hists4Calo_eta[2]->Fill(maxCaloAnaEta);
      }
      if(HLT_AK4CaloJet50_v2){
	hists4Calo_pt[3]->Fill(maxCaloAnaPt);
	hists4Calo_eta[3]->Fill(maxCaloAnaEta);
      }
      if(HLT_AK4CaloJet80_v2){
	hists4Calo_pt[4]->Fill(maxCaloAnaPt);
	hists4Calo_eta[4]->Fill(maxCaloAnaEta);
      }
      if(HLT_AK4CaloJet100_v2){
	hists4Calo_pt[5]->Fill(maxCaloAnaPt);
	hists4Calo_eta[5]->Fill(maxCaloAnaEta);
      }    
    }  
    else{
      if(HLT_AK4CaloJet30_v2) nFake4CaloJet30++;
      if(HLT_AK4CaloJet40_v2) nFake4CaloJet40++;
      if(HLT_AK4CaloJet50_v2) nFake4CaloJet50++;
      if(HLT_AK4CaloJet80_v2) nFake4CaloJet80++;
      if(HLT_AK4CaloJet100_v2) nFake4CaloJet100++;
    }
    

    if(maxPFAnaPt != -1){

      hists4PF_pt[0]->Fill(maxPFAnaPt);
      hists4PF_eta[0]->Fill(maxPFAnaEta);

      if(HLT_AK4PFJet30_v2){
	hists4PF_pt[1]->Fill(maxPFAnaPt);
	hists4PF_eta[1]->Fill(maxPFAnaEta);
      }
      if(HLT_AK4PFJet50_v2){
	hists4PF_pt[2]->Fill(maxPFAnaPt);
	hists4PF_eta[2]->Fill(maxPFAnaEta);
      }
      if(HLT_AK4PFJet80_v2){
	hists4PF_pt[3]->Fill(maxPFAnaPt);
	hists4PF_eta[3]->Fill(maxPFAnaEta);
      }
      if(HLT_AK4PFJet100_v2){
	hists4PF_pt[4]->Fill(maxPFAnaPt);
	hists4PF_eta[4]->Fill(maxPFAnaEta);
      }
    }  
    else{
      if(HLT_AK4PFJet30_v2) nFake4PFJet30++;
      if(HLT_AK4PFJet50_v2) nFake4PFJet50++;
      if(HLT_AK4PFJet80_v2) nFake4PFJet80++;
      if(HLT_AK4PFJet100_v2) nFake4PFJet100++;
    }
  }

  std::cout << "Events matched: " << matched << std::endl;

  std::cout << "True trigger rates: " << std::endl;

  std::cout << "HLT_AK4CaloJet30_v2: " << n4CaloJet30 << ", " << nFake4CaloJet30 << std::endl;
  std::cout << "HLT_AK4CaloJet40_v2: " << n4CaloJet40 << ", " << nFake4CaloJet40 << std::endl;
  std::cout << "HLT_AK4CaloJet50_v2: " << n4CaloJet50 << ", " << nFake4CaloJet50 << std::endl;
  std::cout << "HLT_AK4CaloJet80_v2: " << n4CaloJet80 << ", " << nFake4CaloJet80 << std::endl;
  std::cout << "HLT_AK4CaloJet100_v2: " << n4CaloJet100 << ", " << nFake4CaloJet100 << std::endl;

  std::cout << "HLT_AK4PFJet30_v2: " << n4PFJet30 << ", " << nFake4PFJet30 << std::endl;
  std::cout << "HLT_AK4PFJet50_v2: " << n4PFJet50 << ", " << nFake4PFJet50 << std::endl;
  std::cout << "HLT_AK4PFJet80_v2: " << n4PFJet80 << ", " << nFake4PFJet80 << std::endl;
  std::cout << "HLT_AK4PFJet100_v2: " << n4PFJet100 << ", " << nFake4PFJet100 << std::endl;

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
  TGraphAsymmErrors *a4Calo_pt[nPtTrigCalo], *a4Calo_eta[nPtTrigCalo];
  TGraphAsymmErrors *a4PF_pt[nPtTrigPF], *a4PF_eta[nPtTrigPF];

  for(int i = 0; i < nPtTrigCalo; ++i){
    a4Calo_pt[i] = new TGraphAsymmErrors();
    a4Calo_pt[i]->BayesDivide(hists4Calo_pt[i+1],hists4Calo_pt[0]);
    a4Calo_pt[i]->SetName(ak4CaloName[i]+"_asymm");
    a4Calo_eta[i] = new TGraphAsymmErrors();
    a4Calo_eta[i]->BayesDivide(hists4Calo_eta[i+1],hists4Calo_eta[0]);
    a4Calo_eta[i]->SetName(ak4CaloName[i]+"_eta_asymm");
  }

  for(int i = 0; i < nPtTrigPF; ++i){
    a4PF_pt[i] = new TGraphAsymmErrors();
    a4PF_pt[i]->BayesDivide(hists4PF_pt[i+1],hists4PF_pt[0]);
    a4PF_pt[i]->SetName(ak4PFName[i]+"_asymm");
    a4PF_eta[i] = new TGraphAsymmErrors();
    a4PF_eta[i]->BayesDivide(hists4PF_eta[i+1],hists4PF_eta[0]);
    a4PF_eta[i]->SetName(ak4PFName[i]+"_eta_asymm");
  }

  //save output
  TFile *outFile = TFile::Open("jetTurnOn.root","RECREATE");
  outFile->cd();

  hists4Calo_pt[0]->Write("", TObject::kOverwrite);
  hists4Calo_eta[0]->Write("", TObject::kOverwrite);

  hists4PF_pt[0]->Write("", TObject::kOverwrite);
  hists4PF_eta[0]->Write("", TObject::kOverwrite);

  for(int i = 0; i < nPtTrigCalo; ++i){
    hists4Calo_pt[i+1]->Write("", TObject::kOverwrite);
    hists4Calo_eta[i+1]->Write("", TObject::kOverwrite);
    a4Calo_pt[i]->Write("", TObject::kOverwrite);
    a4Calo_eta[i]->Write("", TObject::kOverwrite);
  }

  for(int i = 0; i < nPtTrigPF; ++i){
    hists4PF_pt[i+1]->Write("", TObject::kOverwrite);
    hists4PF_eta[i+1]->Write("", TObject::kOverwrite);
    a4PF_pt[i]->Write("", TObject::kOverwrite);
    a4PF_eta[i]->Write("", TObject::kOverwrite);
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
