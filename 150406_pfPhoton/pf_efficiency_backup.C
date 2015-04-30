// To validate the PF algorithm this macro the gen tree to the pf tree and counting matched particles and not-matched particles, respectively. The ratio between matched and non-matched particles become efficiency.
// Author : Yeonju Go

//basic c++ header, string ...
#include "../HiForestAnalysis/hiForest.h"
#include "../gammaJetAnalysis/CutAndBinCollection2012.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>
//tree, hist, vector ...
#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>
//canvas, legend, latex ... //cosmetic
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>
//random
#include <TRandom.h>
#include <TStopwatch.h>
#include <ctime>

void pf_efficiency_backup(){
	TString finName="/u/user/goyeonju/PRODUCTION/CMSSW_7_4_0_pre6/src/PFPhoton_hoe025_Forest/merge/ZEE_Hydjet_740pre6_Forest_hoe025_06Apr2015_HiForest.root";
	//TString finName="/u/user/goyeonju/PRODUCTION/CMSSW_7_4_0_pre6/src/PFPhoton_default_Forest/merge/ZEE_Hydjet_740pre6_Forest_hoe050_03Apr2015_HiForest.root";
	const int pfid = 2; // unknown : 0, charged hadron : 1, electron : 2, muon : 3, photon : 4, neutral hadron : 5, h_HF : 6, egamma_HF : 7
	const int pdgid = 11; // electron (11 or -11) , photon 22

	TFile* f1 = new TFile(finName.Data());
	TTree* tgen = (TTree*) f1 -> Get("HiGenParticleAna/hi");
	TTree* tpf = (TTree*) f1 -> Get("pfcandAnalyzer/pfTree");

	const int size = 100000;
	Int_t mult; //gen tree
	Int_t pdg[size];
	Int_t chg[size];
	Float_t pt[size];
	Float_t eta[size];
	Float_t phi[size];
	Int_t nPFpart; //pf tree
	Int_t pfId[size];
	Float_t pfPt[size];
	Float_t pfEta[size];
	Float_t pfPhi[size];	
	TBranch *b_mult; //branch
	TBranch *b_pdg;
	TBranch *b_chg;
	TBranch *b_pt;
	TBranch *b_eta;
	TBranch *b_phi;
	TBranch *b_nPFpart;
	TBranch *b_pfId;
	TBranch *b_pfPt;
	TBranch *b_pfEta;
	TBranch *b_pfPhi;

	tgen -> SetBranchAddress("mult", &mult, &b_mult);
	tgen -> SetBranchAddress("pt", pt, &b_pt);
	tgen -> SetBranchAddress("eta", eta, &b_eta);
	tgen -> SetBranchAddress("phi", phi, &b_phi);
	tgen -> SetBranchAddress("pdg", pdg, &b_pdg);
	tgen -> SetBranchAddress("chg", chg, &b_chg);

	tpf -> SetBranchAddress("nPFpart", &nPFpart, &b_nPFpart);
	tpf -> SetBranchAddress("pfId", pfId, &b_pfId);
	tpf -> SetBranchAddress("pfPt", pfPt, &b_pfPt);
	tpf -> SetBranchAddress("pfEta", pfEta, &b_pfEta);
	tpf -> SetBranchAddress("pfPhi", pfPhi, &b_pfPhi);

	TH1D* h1D_gen = new TH1D("h1D_gen", ";p_{T}^{e};entries", 20,0,200);
	TH1D* h1D_pf  = (TH1D*) h1D_gen->Clone("h1D_pf");
	TH1D* h1D_eff = (TH1D*) h1D_gen->Clone("h1D_eff");
	h1D_eff->GetYaxis()->SetTitle("Efficiency");

	//int Nevt =2;
	int Nevt = tgen->GetEntries();
	for(int ievt=0; ievt< Nevt; ievt++){
		if(ievt%100==0) cout << ">>>>> EVENT " << ievt << " / " << tgen->GetEntries() <<  endl;
		tgen->GetEntry(ievt);
		tpf->GetEntry(ievt);
		for(int igen=0; igen<mult; igen++){
			if(abs(pdg[igen])!=pdgid) continue;
			if(pt[igen]<10) continue;

			float matchedPfPt = 999;
			bool isMatched = 0;
			for(int ipf=0; ipf<nPFpart; ipf++){
				float dR = getDR(eta[igen],phi[igen],pfEta[ipf],pfPhi[ipf]);
				if(dR>0.1) continue;
				if(pfId[ipf]!=pfid) continue;
				isMatched=1;				
				//cout << "what is dR : " << dR << " and pfid : " << pfId[ipf] << endl;
				if(isMatched==1 && abs(pfPt[ipf]-pt[igen])<abs(matchedPfPt-pt[igen])) matchedPfPt = pfPt[ipf];
			}//pf loop
		//	if(ievt < 10) cout << "pt, pfPt : " << pt[igen] << ", " << matchedPfPt<< ", isMatched = " << isMatched <<endl;
			if(isMatched==1) cout << "pt, pfPt : " << pt[igen] << ", " << matchedPfPt <<endl;
			h1D_gen->Fill(pt[igen]);
			if(isMatched==1) h1D_pf->Fill(matchedPfPt);
		}//gen loop
	}//evt loop
	h1D_eff->Divide(h1D_pf, h1D_gen, 1, 1, "B");


	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	///// save the data!
	TCanvas* c1 = new TCanvas("c1", "",400,300);
	c1->cd();
	h1D_eff->Draw();

}
