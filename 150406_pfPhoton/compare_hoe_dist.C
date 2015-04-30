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
#include <TEfficiency.h>
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




void pf_efficiency(TString fname="" ,TH1D* h1D_eff=0, TEfficiency* eff=0, const int pfid=2, const int pdgid=11, float ptThr = 10, const int centLow = 0, const int centUp=10, const float etaLow = 0, const float etaUp = 1.44);
void compare_hoe_dist(const int pfid=2, const int pdgid=11, const int ptThr=10, const int centLow = 0, const int centUp=20, const float etaLow = 0, const float etaUp = 1.44){
	// pfid ;;;;; unknown : 0, charged hadron : 1, electron : 2, muon : 3, photon : 4, neutral hadron : 5, h_HF : 6, egamma_HF : 7
	// pdgid ;;;; electron (11 or -11) , photon 22

	TString fname_default="/u/user/goyeonju/PRODUCTION/CMSSW_7_4_0_pre6/src/PFPhoton_default_Forest/merge/ZEE_Hydjet_740pre6_Forest_hoe050_03Apr2015_HiForest.root";
	TString fname_hoe025="/u/user/goyeonju/PRODUCTION/CMSSW_7_4_0_pre6/src/PFPhoton_hoe025_Forest/merge/ZEE_Hydjet_740pre6_Forest_hoe025_06Apr2015_HiForest.root";

	TH1D* h_default = new TH1D("h_default", ";gen p_{T}^{e};Efficiency", 20,0,200);
	TH1D* h_hoe025= (TH1D*) h_default -> Clone("h_hoe025");
	TEfficiency* eff_default = new TEfficiency("eff_default",";gen p_{T}^{e};Efficiency",20,0,200);
	TEfficiency* eff_hoe025 = new TEfficiency("eff_hoe025",";gen p_{T}^{e};Efficiency",20,0,200);

	pf_efficiency(fname_default, h_default, eff_default, pfid, pdgid, ptThr, centLow, centUp, etaLow, etaUp);
	pf_efficiency(fname_hoe025, h_hoe025, eff_hoe025, pfid, pdgid, ptThr, centLow, centUp, etaLow, etaUp);

	h_default -> SetMarkerColor(1); 
	h_default -> SetAxisRange(0.0,1.2,"Y");
	h_hoe025 -> SetMarkerColor(2); 
	eff_default -> SetMarkerColor(1); 
	eff_hoe025 -> SetMarkerColor(2); 

	TLegend* l1 = new TLegend(0.5, 0.2, 0.8, 0.4);
	l1 -> AddEntry(h_default, "h/e 0.5 (default)");
	l1 -> AddEntry(h_hoe025, "h/e 0.25");

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	///// Draw the eff!
	TCanvas* c1 = new TCanvas("c1", "",300,200);
	c1->cd();
	h_default->Draw("pe");
	h_hoe025->Draw("same pe");
	//	l1->Draw("same");	

	TCanvas* c2 = new TCanvas("c2", "",300,200);
	c2->cd();
	eff_default->Draw("AP");
	eff_hoe025->Draw("AP same");
	
	TString st_id="";
	if(pfid==2) st_id = "electron";
	else if(pfid==4) st_id = "photon";

	TString st_eta="";
	if(etaLow==0) st_eta = "barrel";
	else if(etaLow==1.44) st_eta = "endcap";

	TString st_cent="";
	if(centLow==0) st_cent = "0to10";
	else if(centLow==20) st_cent = "10to30";
	else if(centLow==60) st_cent = "30to100";


	c1->SaveAs(Form("pdf/compare_pf_eff_%s_%s_cent%s.pdf", st_id.Data(), st_eta.Data(), st_cent.Data()));
	c2->SaveAs(Form("pdf/compare_pf_eff_useTEff_%s_%s_cent%s.pdf", st_id.Data(), st_eta.Data(), st_cent.Data()));
} // main function


void pf_efficiency(TString fname, TH1D* h1D_eff, TEfficiency* eff, const int pfid, const int pdgid, const float ptThr, const int centLow, const int centUp, const float etaLow, const float etaUp){
	TFile* f1 = new TFile(fname.Data());
	cout << " Open the file :: " << fname.Data() << endl;
	TTree* tgen = (TTree*) f1 -> Get("HiGenParticleAna/hi");
	TTree* tpf = (TTree*) f1 -> Get("pfcandAnalyzer/pfTree");
	TTree* tevt = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
	tgen->AddFriend(tevt);


	const int size = 100000;
	Int_t hiBin; //evt tree
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
	TBranch *b_hiBin; //branch
	TBranch *b_mult; 
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
	tevt -> SetBranchAddress("hiBin", &hiBin, &b_hiBin);

	TH1D* h1D_pf  = (TH1D*) h1D_eff->Clone("h1D_pf");
	TH1D* h1D_gen = (TH1D*) h1D_eff->Clone("h1D_eff");
	h1D_pf->Sumw2();
	h1D_gen->Sumw2();
	h1D_eff->Sumw2();

	//int Nevt =2;
	int Nevt = tgen->GetEntries();
	for(int ievt=0; ievt< Nevt; ievt++){
		if(ievt%1000==0) cout << ">>>>> EVENT " << ievt << " / " << tgen->GetEntries() <<  endl;
		tgen->GetEntry(ievt);
		tpf->GetEntry(ievt);

		if(!(hiBin>=centLow && hiBin<centUp)) continue; //centrality (0-10), (10-30), (30-100) // hiBin (0,20), (20,60), (60,199);

		for(int igen=0; igen<mult; igen++){
			if(abs(pdg[igen])!=pdgid) continue;
			if(pt[igen]<ptThr) continue;
			if( !(abs(eta[igen])>=etaLow && abs(eta[igen])<etaUp) ) continue; //abs(eta) [0,1.44), [1.44,4) // barrel, endcap

			float matchedPfPt = 999;
			bool isMatched = 0;
			for(int ipf=0; ipf<nPFpart; ipf++){
				float dR = getDR(eta[igen],phi[igen],pfEta[ipf],pfPhi[ipf]);
				if(dR>0.1) continue;
				if(pfId[ipf]!=pfid) continue;
				isMatched=1;
				if(isMatched==1 && abs(pfPt[ipf]-pt[igen])<abs(matchedPfPt-pt[igen])) matchedPfPt = pfPt[ipf];
			}//pf loop
			h1D_gen->Fill(pt[igen]);
			if(isMatched==1) h1D_pf->Fill(pt[igen]);
			//if(isMatched==1) h1D_pf->Fill(matchedPfPt); // this ruins the efficiency. 
			eff->Fill(isMatched,pt[igen]);
		}//gen loop
	}//evt loop
	h1D_eff->Divide(h1D_pf, h1D_gen, 1, 1, "B");
	h1D_eff->SetMarkerStyle(20);
	h1D_eff->SetMarkerSize(1);
}
