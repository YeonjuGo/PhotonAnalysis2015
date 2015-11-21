#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"

#include "EventMatchingCMS.h"
#include "etaPhiFunc.h"

#include <iostream>
#include <string>

const Int_t nRefMax = 500;
const Int_t nGenMax = 100000;

void getLogBins(const Float_t lower, const Float_t higher, const Int_t nBins, Float_t bins[])
{
  Float_t logBins[nBins+1];
  bins[0] = lower;
  bins[nBins] = higher;

  logBins[0] = TMath::Log10(lower);
  logBins[nBins] = TMath::Log10(higher);

  Float_t interval = (logBins[nBins] - logBins[0])/nBins;

  for(Int_t iter = 1; iter < nBins; iter++){
    logBins[iter] = logBins[0] + iter*interval;
    bins[iter] = TMath::Power(10, logBins[iter]);
  }

  return;
}


void FitGauss(TH1F* inHist_p, Float_t &mean, Float_t &meanError, Float_t &res, Float_t &resError)
{
  inHist_p->Fit("gaus", "Q L M", "");

  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);

  Float_t prob = inHist_p->GetFunction("gaus")->GetProb();

  //  if(TMath::Abs(1.00 - mean) < .01) return;                                                                  

  Int_t meanBin = inHist_p->FindBin(mean);
  Float_t meanRMS = -1;
  Float_t total = inHist_p->Integral();

  for(Int_t iter = 0; iter < inHist_p->GetNbinsX(); iter++){
    Int_t lowBound = 0;
    if(meanBin - iter > 0) lowBound = meanBin - iter;

    if(inHist_p->Integral(lowBound, meanBin + iter)/total > .95 || lowBound == 0 || inHist_p->GetBinContent(lowBound) < .01){
      meanRMS = inHist_p->GetBinCenter(meanBin + iter) - inHist_p->GetBinCenter(meanBin);
      break;
    }
  }

  double minPt = inHist_p->GetBinCenter(meanBin) - meanRMS;
  double maxPt = inHist_p->GetBinCenter(meanBin) + meanRMS;

  minPt = std::max(std::max(minPt, 0.0), inHist_p->GetXaxis()->GetXmin());
  maxPt = std::min(maxPt, inHist_p->GetXaxis()->GetXmax());

  inHist_p->Fit("gaus", "Q L M", "", minPt, maxPt);

  if(TMath::Abs(1.00 - mean) < TMath::Abs(1.00 - inHist_p->GetFunction("gaus")->GetParameter(1)) && prob > 0.0001){
    inHist_p->Fit("gaus", "Q L M", "");
    return;
  }

  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);

  return;
}

void basicJetHists(const std::string inName, const std::string outName, const std::string matchName = "")
{
  TH1::SetDefaultSumw2();

  TFile* inFile_p = new TFile(inName.c_str(), "READ");
  TTree* jet4Tree_p = (TTree*)inFile_p->Get("akPu4CaloJetAnalyzer/t");
  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  TTree* genTree_p = (TTree*)inFile_p->Get("HiGenParticleAna/hi");

  Int_t hi_evt;
  Int_t hi_lumi;

  Int_t nref_;
  Float_t jtpt_[nRefMax];
  Float_t jtphi_[nRefMax];
  Float_t jteta_[nRefMax];

  Float_t refpt_[nRefMax];
  Float_t refphi_[nRefMax];
  Float_t refeta_[nRefMax];

  Int_t ngen_;
  Float_t genpt_[nRefMax];
  Float_t genphi_[nRefMax];
  Float_t geneta_[nRefMax];

  Int_t mult_ = 0;
  Float_t pt_[nGenMax];
  Int_t pdg_[nGenMax];

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("evt", 1);
  hiTree_p->SetBranchStatus("lumi", 1);

  hiTree_p->SetBranchAddress("evt", &hi_evt);
  hiTree_p->SetBranchAddress("lumi", &hi_lumi);

  jet4Tree_p->SetBranchStatus("*", 0);
  jet4Tree_p->SetBranchStatus("nref", 1);
  jet4Tree_p->SetBranchStatus("jtpt", 1);
  jet4Tree_p->SetBranchStatus("jtphi", 1);
  jet4Tree_p->SetBranchStatus("jteta", 1);

  jet4Tree_p->SetBranchStatus("refpt", 1);
  jet4Tree_p->SetBranchStatus("refphi", 1);
  jet4Tree_p->SetBranchStatus("refeta", 1);

  jet4Tree_p->SetBranchStatus("ngen", 1);
  jet4Tree_p->SetBranchStatus("genpt", 1);
  jet4Tree_p->SetBranchStatus("genphi", 1);
  jet4Tree_p->SetBranchStatus("geneta", 1);

  jet4Tree_p->SetBranchAddress("nref", &nref_);
  jet4Tree_p->SetBranchAddress("jtpt", jtpt_);
  jet4Tree_p->SetBranchAddress("jtphi", jtphi_);
  jet4Tree_p->SetBranchAddress("jteta", jteta_);

  jet4Tree_p->SetBranchAddress("refpt", refpt_);
  jet4Tree_p->SetBranchAddress("refphi", refphi_);
  jet4Tree_p->SetBranchAddress("refeta", refeta_);

  jet4Tree_p->SetBranchAddress("ngen", &ngen_);
  jet4Tree_p->SetBranchAddress("genpt", genpt_);
  jet4Tree_p->SetBranchAddress("genphi", genphi_);
  jet4Tree_p->SetBranchAddress("geneta", geneta_);

  
  genTree_p->SetBranchStatus("*", 0);
  /*
  genTree_p->SetBranchStatus("mult", 1);
  genTree_p->SetBranchStatus("pt", 1);
  genTree_p->SetBranchStatus("pdg", 1);

  genTree_p->SetBranchAddress("mult", &mult_);
  genTree_p->SetBranchAddress("pt", &pt_);
  genTree_p->SetBranchAddress("pdg", &pdg_);
  */

  Int_t nMuBins = 5;
  Float_t muBin[nMuBins+1] = {4.0, 6.0, 8.5, 11.0, 13.5, 16.0};
  TH1F* muRateHist_h = new TH1F("muRateHist_h", "muRateHist_h", nMuBins, muBin);
  TH1F* muJetRateHist_h = new TH1F("muJetRateHist_h", "muJetRateHist_h", nMuBins, muBin);
  TH1F* muGammaRateHist_h = new TH1F("muGammaRateHist_h", "muGammaRateHist_h", nMuBins, muBin);

  TFile* matchFile_p;
  TTree* matchTree_p;

  Int_t match_lumi;
  Int_t match_evt;
  Int_t jetPt_[8];
  Int_t jetHwPhi_[8];
  Int_t jetHwEta_[8];

  const Float_t phiMap[18] = {0, .3491, .6981, 1.047, 1.396, 1.745, 2.094, 2.443, 2.793, 3.1415, -2.793, -2.443, -2.094, -1.745, -1.396, -1.047, -.6981, -.3491};
  const Float_t etaMap[22] = {-3.2, -3.1, -3.0, -2.9, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 2.9, 3.0, 3.1, 3.2};

  EventMatchingCMS *matcher = new EventMatchingCMS();

  if(strcmp(matchName.c_str(), "") != 0){
    matchFile_p = new TFile(matchName.c_str(), "READ");
    matchTree_p = (TTree*)matchFile_p->Get("L1UpgradeTree");

    matchTree_p->SetBranchStatus("*", 0);
    matchTree_p->SetBranchStatus("lumi", 1);
    matchTree_p->SetBranchStatus("event", 1);
    matchTree_p->SetBranchStatus("jet_pt", 1);
    matchTree_p->SetBranchStatus("jet_hwPhi", 1);
    matchTree_p->SetBranchStatus("jet_hwEta", 1);

    matchTree_p->SetBranchAddress("lumi", &match_lumi);
    matchTree_p->SetBranchAddress("event", &match_evt);
    matchTree_p->SetBranchAddress("jet_pt", jetPt_);
    matchTree_p->SetBranchAddress("jet_hwPhi", jetHwPhi_);
    matchTree_p->SetBranchAddress("jet_hwEta", jetHwEta_);

    for(Long64_t entry = 0; entry < matchTree_p->GetEntries(); ++entry){
      matchTree_p->GetEntry(entry);
      matcher->addEvent(match_evt, match_lumi, 0, entry);
    }
  }


  const std::string genLeadJtName[3] = {"genLeadJet4Pt_h", "genLeadJet4Eta_h", "genLeadJet4Phi_h"};
  const std::string recoLeadJtName[3] = {"recoLeadJet4Pt_h", "recoLeadJet4Eta_h", "recoLeadJet4Phi_h"};
  const Int_t bins[3] = {100, 50, 50};
  const Float_t lower[3] = {-0.5, -5.0, (Float_t)-TMath::Pi()};
  const Float_t upper[3] = {199.5, 5.0, (Float_t)TMath::Pi()};

  TH1F* genLeadJet4_p[3];
  TH1F* recoLeadJet4_p[3];
  
  for(Int_t iter = 0; iter < 3; iter++){
    genLeadJet4_p[iter] = new TH1F(genLeadJtName[iter].c_str(), genLeadJtName[iter].c_str(), bins[iter], lower[iter], upper[iter]);
    recoLeadJet4_p[iter] = new TH1F(recoLeadJtName[iter].c_str(), recoLeadJtName[iter].c_str(), bins[iter], lower[iter], upper[iter]);
  }

  const Int_t nJECBin = 10;
  Float_t jecPtBins[nJECBin+1];
  getLogBins(20, 200, nJECBin, jecPtBins);
  Float_t jecEtaBins[nJECBin+1];
  for(Int_t iter = 0; iter < nJECBin+1; iter++){
    jecEtaBins[iter] = (-5.0 + iter*10.0/nJECBin);
    std::cout << jecEtaBins[iter] << std::endl;
  }

  std::cout << "JEC BINS" << std::endl;
  for(Int_t iter = 0; iter < nJECBin + 1; iter++){
    std::cout << jecPtBins[iter] << std::endl;
  }

  TH1F* leadJetJEC_Pt_h = new TH1F("leadJetJEC_Pt_h", "leadJetJEC_Pt_h", nJECBin, jecPtBins);
  TH1F* leadJetJEC_Eta_h = new TH1F("leadJetJEC_Eta_h", "leadJetJEC_Eta_h", nJECBin, jecEtaBins);
  TH1F* leadJetFake_h = new TH1F("leadJetFake_h", "leadJetFake_h", nJECBin, jecPtBins);
  TH1F* leadJetFake2_h = new TH1F("leadJetFake2_h", "leadJetFake2_h", nJECBin, jecPtBins);
  TH1F* leadJetFake3_h = new TH1F("leadJetFake3_h", "leadJetFake3_h", nJECBin, jecPtBins);
  TH1F* leadJetL1Scale_h = new TH1F("leadJetL1Scale_h", "leadJetL1Scale_h", nJECBin, jecPtBins);


  TH1F* leadJetJEC_Pt_Bin_h[nJECBin];
  TH1F* leadJetJEC_Eta_Bin_h[nJECBin];
  TH1F* leadJetFake_Bin_h[nJECBin];
  TH1F* leadJetFake2_Bin_h[nJECBin];
  TH1F* leadJetFake3_Bin_h[nJECBin];
  TH1F* leadJetL1Scale_Bin_h[nJECBin];
  for(Int_t iter = 0; iter < nJECBin; iter++){
    leadJetJEC_Pt_Bin_h[iter] = new TH1F(Form("leadJetJEC_Pt_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), Form("leadJetJEC_Pt_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), 100, 0, 4);
    leadJetJEC_Eta_Bin_h[iter] = new TH1F(Form("leadJetJEC_Eta_%d_h", iter), Form("leadJetJEC_Eta_%d_h", iter), 100, 0, 4);
    leadJetFake_Bin_h[iter] = new TH1F(Form("leadJetFake_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), Form("leadJetFake_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), 2, -0.5, 1.5);
    leadJetFake2_Bin_h[iter] = new TH1F(Form("leadJetFake2_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), Form("leadJetFake2_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), 2, -0.5, 1.5);
    leadJetFake3_Bin_h[iter] = new TH1F(Form("leadJetFake3_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), Form("leadJetFake3_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), 2, -0.5, 1.5);
    leadJetL1Scale_Bin_h[iter] = new TH1F(Form("leadJetL1Scale_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), Form("leadJetL1Scale_%d_%d_h", (Int_t)jecPtBins[iter], (Int_t)jecPtBins[iter+1]), 100, 0, 4);
  }


  Int_t nEntries = jet4Tree_p->GetEntries();
  Int_t matched = 0;
  Int_t superFake = 0;

  std::cout << "Entries: " << nEntries << std::endl;

  std::cout << "SCALING HISTS" << std::endl;

  if(strcmp(matchName.c_str(), "") != 0){
    for(Int_t entry = 0; entry < nEntries; entry++){
      if(entry%10000 == 0) std::cout << "Event: " << entry << std::endl;
      jet4Tree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);

      if(strcmp(matchName.c_str(), "") != 0){
	long long hlt_entry = matcher->retrieveEvent(hi_evt, hi_lumi, 0);
	if(hlt_entry == -1) continue;
	matchTree_p->GetEntry(hlt_entry);
	matched++;
	
	if(jetPt_[0] > 0){
	  Float_t tempPhi = phiMap[jetHwPhi_[0]];
	  Float_t tempEta = etaMap[jetHwEta_[0]];
	  Int_t tempJtPos = -1;
	  
	  for(Int_t iter = 0; iter < nref_; iter++){
	    if(getDR(tempEta, tempPhi, jteta_[iter], jtphi_[iter]) < 1.0){
	      tempJtPos = iter;
	      break;
	    }
	  }
	  
	  if(tempJtPos != -1)
	    if(refpt_[tempJtPos] <= 0) tempJtPos = -1;	

	  if(tempJtPos != -1){	  
	    for(Int_t iter = 0; iter < nJECBin; iter++){
	      if(jetPt_[0] < jecPtBins[iter+1]){
		leadJetL1Scale_Bin_h[iter]->Fill(jtpt_[tempJtPos]/jetPt_[0]);
		break;
	      }
	    }
	  }

	}
      }
      
    }
  }


  delete matcher;
  matcher = new EventMatchingCMS();

  if(strcmp(matchName.c_str(), "") != 0){
    matchFile_p = new TFile(matchName.c_str(), "READ");
    matchTree_p = (TTree*)matchFile_p->Get("L1UpgradeTree");

    matchTree_p->SetBranchStatus("*", 0);
    matchTree_p->SetBranchStatus("lumi", 1);
    matchTree_p->SetBranchStatus("event", 1);
    matchTree_p->SetBranchStatus("jet_pt", 1);
    matchTree_p->SetBranchStatus("jet_hwPhi", 1);
    matchTree_p->SetBranchStatus("jet_hwEta", 1);

    matchTree_p->SetBranchAddress("lumi", &match_lumi);
    matchTree_p->SetBranchAddress("event", &match_evt);
    matchTree_p->SetBranchAddress("jet_pt", jetPt_);
    matchTree_p->SetBranchAddress("jet_hwPhi", jetHwPhi_);
    matchTree_p->SetBranchAddress("jet_hwEta", jetHwEta_);

    for(Long64_t entry = 0; entry < matchTree_p->GetEntries(); ++entry){
      matchTree_p->GetEntry(entry);
      matcher->addEvent(match_evt, match_lumi, 0, entry);
    }
  }

  Int_t nMu5 = 0;
  Int_t nMu7 = 0;
  Int_t nMu10 = 0;
  Int_t nMu12 = 0;
  Int_t nMu15 = 0;

  Int_t nMu5_Jet = 0;
  Int_t nMu7_Jet = 0;
  Int_t nMu10_Jet = 0;
  Int_t nMu12_Jet = 0;
  Int_t nMu15_Jet = 0;

  Int_t nMu5_Gamma = 0;
  Int_t nMu7_Gamma = 0;
  Int_t nMu10_Gamma = 0;
  Int_t nMu12_Gamma = 0;
  Int_t nMu15_Gamma = 0;

  for(Int_t entry = 0; entry < nEntries; entry++){
    if(entry%10000 == 0) std::cout << "Event: " << entry << std::endl;
    jet4Tree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);
    genTree_p->GetEntry(entry);

    Bool_t mu5 = false;
    Bool_t mu7 = false;
    Bool_t mu10 = false;
    Bool_t mu12 = false;
    Bool_t mu15 = false;

    Bool_t gamma15 = false;

    for(Int_t iter = 0; iter < mult_; iter++){
      if(TMath::Abs(pdg_[iter]) == 13 && !mu15){
	if(pt_[iter] > 15) mu15 = true;
	else if(pt_[iter] > 12) mu12 = true;
	else if(pt_[iter] > 10) mu10 = true;
	else if(pt_[iter] > 7) mu7 = true;
	else if(pt_[iter] > 5) mu5 = true;	  
      }
      else if(TMath::Abs(pdg_[iter]) == 22){
	if(pt_[iter] > 15) gamma15 = true;
      }

      if(mu15 && gamma15) break;
    }

    if(mu15){
      nMu5++;
      nMu7++;
      nMu10++;
      nMu12++;
      nMu15++;

      if(gamma15){
	nMu5_Gamma++;
	nMu7_Gamma++;
	nMu10_Gamma++;
	nMu12_Gamma++;
	nMu15_Gamma++;
      }
    }
    else if(mu12){
      nMu5++;
      nMu7++;
      nMu10++;
      nMu12++;

      if(gamma15){
        nMu5_Gamma++;
        nMu7_Gamma++;
        nMu10_Gamma++;
        nMu12_Gamma++;
      }
    }
    else if(mu10){
      nMu5++;
      nMu7++;
      nMu10++;

      if(gamma15){
        nMu5_Gamma++;
        nMu7_Gamma++;
        nMu10_Gamma++;
      }
    }
    else if(mu7){
      nMu5++;
      nMu7++;

      if(gamma15){
        nMu5_Gamma++;
        nMu7_Gamma++;
      }
    }
    else if(mu5){
      nMu5++;
      if(gamma15) nMu5_Gamma++;
    }

    if(nref_ != 0){
      recoLeadJet4_p[0]->Fill(jtpt_[0]);
      recoLeadJet4_p[1]->Fill(jteta_[0]);
      recoLeadJet4_p[2]->Fill(jtphi_[0]);

      if(refpt_[0] > 0 && TMath::Abs(refeta_[0]) < 5.0){
	for(Int_t iter = 0; iter < nJECBin; iter++){
	  if(refpt_[0] < jecPtBins[iter+1]){
	    leadJetJEC_Pt_Bin_h[iter]->Fill(jtpt_[0]/refpt_[0]);
	    break;
	  }
	}

	if(refpt_[0] > 50){
	  for(Int_t iter = 0; iter < nJECBin; iter++){
	    if(refeta_[0] < jecEtaBins[iter+1]){
	      leadJetJEC_Eta_Bin_h[iter]->Fill(jtpt_[0]/refpt_[0]);
	      break;
	    }
	  }
	}
      }

      if(TMath::Abs(jteta_[0]) < 2.0){
	for(Int_t iter = 0; iter < nJECBin; iter++){
	  if(jtpt_[0] < jecPtBins[iter+1]){
	    if(refpt_[0] > 0) leadJetFake_Bin_h[iter]->Fill(0);
	    else leadJetFake_Bin_h[iter]->Fill(1);
	    break;
	  }
	}
      }

      if(mu15 && jtpt_[0] > 60){
	nMu5_Jet++;
	nMu7_Jet++;
	nMu10_Jet++;
	nMu12_Jet++;
	nMu15_Jet++;
      }
      else if(mu12 && jtpt_[0] > 60){
	nMu5_Jet++;
	nMu7_Jet++;
	nMu10_Jet++;
	nMu12_Jet++;
      }
      else if(mu10 && jtpt_[0] > 60){
	nMu5_Jet++;
	nMu7_Jet++;
	nMu10_Jet++;
      }
      else if(mu7 && jtpt_[0] > 60){
	nMu5_Jet++;
	nMu7_Jet++;
      }
      else if(mu5 && jtpt_[0] > 60) nMu5_Jet++;

    }

    if(ngen_ != 0){
      Float_t maxGenPt = -1;
      Float_t maxGenEta = -1;
      Float_t maxGenPhi = -1;

      for(Int_t iter = 0; iter < ngen_; iter++){
	if(genpt_[iter] > maxGenPt){
	  maxGenPt = genpt_[iter];
	  maxGenPhi = genphi_[iter];
	  maxGenEta = geneta_[iter];
	}
      }

      genLeadJet4_p[0]->Fill(maxGenPt);
      genLeadJet4_p[1]->Fill(maxGenEta);
      genLeadJet4_p[2]->Fill(maxGenPhi);
    }

    if(strcmp(matchName.c_str(), "") != 0){
      long long hlt_entry = matcher->retrieveEvent(hi_evt, hi_lumi, 0);
      if(hlt_entry == -1) continue;
      matchTree_p->GetEntry(hlt_entry);
      matched++;

      if(jetPt_[0] > 0){
	Float_t tempPhi = phiMap[jetHwPhi_[0]];
	Float_t tempEta = etaMap[jetHwEta_[0]];
	Int_t tempJtPos = -1;

	for(Int_t iter = 0; iter < nref_; iter++){
	  if(getDR(tempEta, tempPhi, jteta_[iter], jtphi_[iter]) < 0.8){
	    tempJtPos = iter;
	    break;
	  }
	}
	
	if(tempJtPos != -1){
	  Int_t tempJtPos2 = tempJtPos;
	  if(refpt_[tempJtPos] <= 0) tempJtPos = -1;
	  
 
	  for(Int_t iter = 0; iter < nJECBin; iter++){
	    if(jtpt_[tempJtPos2] < jecPtBins[iter+1]){
	      if(tempJtPos == -1) leadJetFake3_Bin_h[iter]->Fill(1);
	      else leadJetFake3_Bin_h[iter]->Fill(0);
	      break;
	    }
	  }
	}
	else{
	  superFake++;
	}

	for(Int_t iter = 0; iter < nJECBin; iter++){
	  if(jetPt_[0] < jecPtBins[iter+1]){
	    jetPt_[0] *= leadJetL1Scale_Bin_h[iter]->GetMean(); 
	    break;
	  }
	}

	
	for(Int_t iter = 0; iter < nJECBin; iter++){
          if(jetPt_[0] < jecPtBins[iter+1]){
            if(tempJtPos == -1) leadJetFake2_Bin_h[iter]->Fill(1);
	    else leadJetFake2_Bin_h[iter]->Fill(0);
            break;
          }
	}
      }
    }
  }

  std::cout << "Events Matched: " << matched << std::endl;
  std::cout << "superFake: " << superFake << std::endl;

  std::cout << "Mu15: " << nMu15 << std::endl;
  std::cout << "Mu12: " << nMu12 << std::endl;
  std::cout << "Mu10: " << nMu10 << std::endl;
  std::cout << "Mu7: " << nMu7 << std::endl;
  std::cout << "Mu5: " << nMu5 << std::endl;

  muRateHist_h->SetBinContent(1, nMu5);
  muRateHist_h->SetBinContent(2, nMu7);
  muRateHist_h->SetBinContent(3, nMu10);
  muRateHist_h->SetBinContent(4, nMu12);
  muRateHist_h->SetBinContent(5, nMu15);

  muRateHist_h->SetBinError(1, TMath::Sqrt(nMu5));
  muRateHist_h->SetBinError(2, TMath::Sqrt(nMu7));
  muRateHist_h->SetBinError(3, TMath::Sqrt(nMu10));
  muRateHist_h->SetBinError(4, TMath::Sqrt(nMu12));
  muRateHist_h->SetBinError(5, TMath::Sqrt(nMu15));

  muRateHist_h->Scale(30000.0/nEntries);

  std::cout << "MuGamma15: " << nMu15_Gamma << std::endl;
  std::cout << "MuGamma12: " << nMu12_Gamma << std::endl;
  std::cout << "MuGamma10: " << nMu10_Gamma << std::endl;
  std::cout << "MuGamma7: " << nMu7_Gamma << std::endl;
  std::cout << "MuGamma5: " << nMu5_Gamma << std::endl;

  muGammaRateHist_h->SetBinContent(1, nMu5_Gamma);
  muGammaRateHist_h->SetBinContent(2, nMu7_Gamma);
  muGammaRateHist_h->SetBinContent(3, nMu10_Gamma);
  muGammaRateHist_h->SetBinContent(4, nMu12_Gamma);
  muGammaRateHist_h->SetBinContent(5, nMu15_Gamma);

  muGammaRateHist_h->SetBinError(1, TMath::Sqrt(nMu5_Gamma));
  muGammaRateHist_h->SetBinError(2, TMath::Sqrt(nMu7_Gamma));
  muGammaRateHist_h->SetBinError(3, TMath::Sqrt(nMu10_Gamma));
  muGammaRateHist_h->SetBinError(4, TMath::Sqrt(nMu12_Gamma));
  muGammaRateHist_h->SetBinError(5, TMath::Sqrt(nMu15_Gamma));

  muGammaRateHist_h->Scale(30000.0/nEntries);

  std::cout << "Mu15_Jet: " << nMu15_Jet << std::endl;
  std::cout << "Mu12_Jet: " << nMu12_Jet << std::endl;
  std::cout << "Mu10_Jet: " << nMu10_Jet << std::endl;
  std::cout << "Mu7_Jet: " << nMu7_Jet << std::endl;
  std::cout << "Mu5_Jet: " << nMu5_Jet << std::endl;

  muJetRateHist_h->SetBinContent(1, nMu5_Jet);
  muJetRateHist_h->SetBinContent(2, nMu7_Jet);
  muJetRateHist_h->SetBinContent(3, nMu10_Jet);
  muJetRateHist_h->SetBinContent(4, nMu12_Jet);
  muJetRateHist_h->SetBinContent(5, nMu15_Jet);

  muJetRateHist_h->SetBinError(1, TMath::Sqrt(nMu5_Jet));
  muJetRateHist_h->SetBinError(2, TMath::Sqrt(nMu7_Jet));
  muJetRateHist_h->SetBinError(3, TMath::Sqrt(nMu10_Jet));
  muJetRateHist_h->SetBinError(4, TMath::Sqrt(nMu12_Jet));
  muJetRateHist_h->SetBinError(5, TMath::Sqrt(nMu15_Jet));

  muJetRateHist_h->Scale(30000.0/nEntries);

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");

  muRateHist_h->Write("", TObject::kOverwrite);  
  muJetRateHist_h->Write("", TObject::kOverwrite);  
  muGammaRateHist_h->Write("", TObject::kOverwrite);  

  for(Int_t iter = 0; iter < 3; iter++){
    recoLeadJet4_p[iter]->Scale(1./nEntries);
    recoLeadJet4_p[iter]->Write("", TObject::kOverwrite);

    genLeadJet4_p[iter]->Scale(1./nEntries);
    genLeadJet4_p[iter]->Write("", TObject::kOverwrite);
  }
  for(Int_t iter = 0; iter < nJECBin; iter++){
    Float_t mean = 0;
    Float_t meanErr = 0;
    Float_t res = 0;
    Float_t resErr = 0;

    FitGauss(leadJetJEC_Pt_Bin_h[iter], mean, meanErr, res, resErr);
    leadJetJEC_Pt_h->SetBinContent(iter+1, mean);
    leadJetJEC_Pt_h->SetBinError(iter+1, meanErr);
    leadJetJEC_Pt_Bin_h[iter]->Write("", TObject::kOverwrite);

    FitGauss(leadJetJEC_Eta_Bin_h[iter], mean, meanErr, res, resErr);
    leadJetJEC_Eta_h->SetBinContent(iter+1, mean);
    leadJetJEC_Eta_h->SetBinError(iter+1, meanErr);
    leadJetJEC_Eta_Bin_h[iter]->Write("", TObject::kOverwrite);

    leadJetFake_h->SetBinContent(iter+1, leadJetFake_Bin_h[iter]->GetMean());
    leadJetFake_h->SetBinError(iter+1, leadJetFake_Bin_h[iter]->GetMeanError());
    leadJetFake_Bin_h[iter]->Write("", TObject::kOverwrite);

    if(strcmp(matchName.c_str(), "") != 0){
      leadJetFake2_h->SetBinContent(iter+1, leadJetFake2_Bin_h[iter]->GetMean());
      leadJetFake2_h->SetBinError(iter+1, leadJetFake2_Bin_h[iter]->GetMeanError());
      leadJetFake2_Bin_h[iter]->Write("", TObject::kOverwrite);

      leadJetFake3_h->SetBinContent(iter+1, leadJetFake3_Bin_h[iter]->GetMean());
      leadJetFake3_h->SetBinError(iter+1, leadJetFake3_Bin_h[iter]->GetMeanError());
      leadJetFake3_Bin_h[iter]->Write("", TObject::kOverwrite);

      leadJetL1Scale_h->SetBinContent(iter+1, leadJetL1Scale_Bin_h[iter]->GetMean());
      leadJetL1Scale_h->SetBinError(iter+1, leadJetL1Scale_Bin_h[iter]->GetMeanError());
      leadJetL1Scale_Bin_h[iter]->Write("", TObject::kOverwrite);

      std::cout << iter << ": " << leadJetL1Scale_Bin_h[iter]->GetMean() << std::endl;
    }
  }
  leadJetJEC_Pt_h->Write("", TObject::kOverwrite);
  leadJetJEC_Eta_h->Write("", TObject::kOverwrite);
  leadJetFake_h->Write("", TObject::kOverwrite);
  if(strcmp(matchName.c_str(), "") != 0){
    leadJetFake2_h->Write("", TObject::kOverwrite);
    leadJetFake3_h->Write("", TObject::kOverwrite);
    leadJetL1Scale_h->Write("", TObject::kOverwrite);
  }

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < 3; iter++){
    delete recoLeadJet4_p[iter];    
    delete genLeadJet4_p[iter];    
  }

  for(Int_t iter = 0; iter < nJECBin; iter++){
    delete leadJetJEC_Pt_Bin_h[iter];
    delete leadJetJEC_Eta_Bin_h[iter];
    delete leadJetFake_Bin_h[iter];

    if(strcmp(matchName.c_str(), "") != 0){
      delete leadJetFake2_Bin_h[iter];
      delete leadJetFake3_Bin_h[iter];
      delete leadJetL1Scale_Bin_h[iter];
    }
  }

  delete leadJetJEC_Pt_h;

  delete leadJetJEC_Eta_h;

  delete leadJetFake_h;
  delete leadJetFake2_h;
  delete leadJetFake3_h;
  delete leadJetL1Scale_h;

  delete inFile_p;

  return;
}
