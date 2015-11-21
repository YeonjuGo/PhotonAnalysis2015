#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <iostream>

void makeGenSkim()
{
  TFile* inFile_p = new TFile("HydjetMB_502TeV_740pre8_MCHI2_74_V3_rctconfigNoCuts_HiForestAndEmulatorAndHLT_v7.root", "READ");
  TTree* genTree_p = (TTree*)inFile_p->Get("HiGenParticleAna/hi");
  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");

  Int_t maxMult_ = 100000;
  Int_t mult_;
  Int_t pdg_[maxMult_], chg_[maxMult_];
  Float_t pt_[maxMult_], phi_[maxMult_], eta_[maxMult_];

  Int_t hiBin_;
  Float_t vx_, vy_, vz_;

  Int_t pcollisionEventSelection_;

  genTree_p->SetBranchStatus("*", 0);
  genTree_p->SetBranchStatus("mult", 1);
  genTree_p->SetBranchStatus("pdg", 1);
  genTree_p->SetBranchStatus("chg", 1);
  genTree_p->SetBranchStatus("pt", 1);
  genTree_p->SetBranchStatus("phi", 1);
  genTree_p->SetBranchStatus("eta", 1);

  genTree_p->SetBranchAddress("mult", &mult_);
  genTree_p->SetBranchAddress("pdg", pdg_);
  genTree_p->SetBranchAddress("chg", chg_);
  genTree_p->SetBranchAddress("pt", pt_);
  genTree_p->SetBranchAddress("phi", phi_);
  genTree_p->SetBranchAddress("eta", eta_);

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("hiBin", 1);
  hiTree_p->SetBranchStatus("vx", 1);
  hiTree_p->SetBranchStatus("vy", 1);
  hiTree_p->SetBranchStatus("vz", 1);

  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  hiTree_p->SetBranchAddress("vx", &vx_);
  hiTree_p->SetBranchAddress("vy", &vy_);
  hiTree_p->SetBranchAddress("vz", &vz_);

  skimTree_p->SetBranchStatus("*", 0);
  skimTree_p->SetBranchStatus("pcollisionEventSelection", 1);

  skimTree_p->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);

  TFile* outFile_p = new TFile("HydjetMB_502_BackgroundSkim.root", "RECREATE");
  TTree* outTree_p = new TTree("genHydTree", "genHydTree");

  Int_t nGen_;
  Int_t genID_[maxMult_], genChg_[maxMult_];
  Float_t genPt_[maxMult_], genPhi_[maxMult_], genEta_[maxMult_];

  outTree_p->Branch("nGen", &nGen_, "nGen/I");
  outTree_p->Branch("genID", genID_, "genID[nGen]/I");
  outTree_p->Branch("genChg", genChg_, "genChg[nGen]/I");
  outTree_p->Branch("genPt", genPt_, "genPt[nGen]/F");
  outTree_p->Branch("genPhi", genPhi_, "genPhi[nGen]/F");
  outTree_p->Branch("genEta", genEta_, "genEta[nGen]/F");
  outTree_p->Branch("hiBin", &hiBin_, "hiBin/I");
  outTree_p->Branch("vx", &vx_, "vx/F");
  outTree_p->Branch("vy", &vy_, "vy/F");
  outTree_p->Branch("vz", &vz_, "vz/F");

  Int_t nEntries = genTree_p->GetEntries();

  for(Int_t iter = 0; iter < nEntries; iter++){
    if(iter%10000 == 0) std::cout << "Entry: " << iter << std::endl;

    genTree_p->GetEntry(iter);
    hiTree_p->GetEntry(iter);
    skimTree_p->GetEntry(iter);

    if(TMath::Abs(vz_) > 15) continue;
    if(!pcollisionEventSelection_) continue;

    nGen_ = 0;

    for(Int_t iter2 = 0; iter2 < mult_; iter2++){
      if(TMath::Abs(eta_[iter2]) > 2.4) continue;
      if(pt_[iter2] < 0.5) continue;
      if(chg_[iter2] == 0) continue;

      genID_[nGen_] = pdg_[iter2];
      genChg_[nGen_] = chg_[iter2];
      genPt_[nGen_] = pt_[iter2];
      genPhi_[nGen_] = phi_[iter2];
      genEta_[nGen_] = eta_[iter2];
      nGen_++;
    }

    outTree_p->Fill();

    if(outTree_p->GetEntries() == 50000) break;
  }

  outFile_p->cd();
  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;
  outFile_p->Close();
  delete outFile_p;

  delete inFile_p;

  return;
}
