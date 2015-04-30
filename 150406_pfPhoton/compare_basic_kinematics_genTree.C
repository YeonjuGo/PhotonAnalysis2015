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

void compareTwo(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, TCut theCut="", int numbering=1);
void compare_basic_kinematics_genTree(){

	const int pfid = 2; // unknown : 0, charged hadron : 1, electron : 2, muon : 3, photon : 4, neutral hadron : 5, h_HF : 6, egamma_HF : 7
	const int pdgid = 11; // electron (11 or -11) , photon 22
	const float ptThr= 10; // electron (11 or -11) , photon 22

	TString fname_default="/u/user/goyeonju/PRODUCTION/CMSSW_7_4_0_pre6/src/PFPhoton_default_Forest/merge/ZEE_Hydjet_740pre6_Forest_hoe050_03Apr2015_HiForest.root";
	TString fname_hoe025="/u/user/goyeonju/PRODUCTION/CMSSW_7_4_0_pre6/src/PFPhoton_hoe025_Forest/merge/ZEE_Hydjet_740pre6_Forest_hoe025_06Apr2015_HiForest.root";
 
        TFile* f1 = new TFile(fname_default.Data());
        TTree* t1 = (TTree*) f1 -> Get("HiGenParticleAna/hi");
        TFile* f2 = new TFile(fname_hoe025.Data());
        TTree* t2 = (TTree*) f2 -> Get("HiGenParticleAna/hi");

	compareTwo(t1, t2, "mult",50,0,60000,"");
	compareTwo(t1, t2, "pdg",60,-30,30,"");
	compareTwo(t1, t2, "pt",50,0,200,"abs(pdg)==11",1);
	compareTwo(t1, t2, "pt",50,0,200,"abs(pdg)==11&&abs(eta)<1.44",2);
	compareTwo(t1, t2, "pt",50,0,200,"abs(pdg)==11&&abs(eta)>1.44",3);
	compareTwo(t1, t2, "pt",50,0,200,"abs(pdg)==22",4);
	compareTwo(t1, t2, "pt",50,0,200,"abs(pdg)==22&&abs(eta)<1.44",5);
	compareTwo(t1, t2, "pt",50,0,200,"abs(pdg)==22&&abs(eta)>1.44",6);
#if 0
	compareTwo(t1, t2, "pfId",8, 0, 8, "",0);
	compareTwo(t1, t2, "pfId",8, 0, 8, "abs(pfEta)<1.44&&pfPt>10",1);
	compareTwo(t1, t2, "pfPt",50, 0, 200, "pfId==2&&abs(pfEta)<1.44",0); //electron barrel
	compareTwo(t1, t2, "pfPt",50, 0, 200, "pfId==2&&abs(pfEta)>1.44",1); //electron endcap
	compareTwo(t1, t2, "pfPt",50, 0, 200, "pfId==4&&abs(pfEta)<1.44",2); //photon barrel
	compareTwo(t1, t2, "pfPt",50, 0, 200, "pfId==4&&abs(pfEta)>1.44",3); //photon endcap
#endif


} // main function

void compareTwo(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, TCut theCut, int numbering)  {
	TCanvas* c=  new TCanvas(Form("c_%s_%d",var.Data(),numbering),"", 400,800);
	c->Divide(1,2);
	c->cd(1);
	gPad->SetLogy();
	TH1D* h1 = new TH1D(Form("h1_%s_%d",var.Data(),numbering), Form(";%s;",var.Data()), nBins,xMin,xMax);
	TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%d",var.Data(),numbering));
	h1->Sumw2();
	h2->Sumw2();
	t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
	t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);	
	h1->Scale( 1. / t1->GetEntries());
	h2->Scale( 1. / t2->GetEntries());
	h2->SetMarkerStyle(20);
	h2->SetMarkerSize(0.8);
	double range = cleverRange(h1,h2,1.5,1.e-4);
	h1->DrawCopy("hist");
	h2->DrawCopy("L same");
	c->cd(2);
	h2->Divide(h1);
	h2->SetYTitle("(h/e 0.25)/(h/e 0.5) Ratio");

	double ratioRange = getCleverRange(h2);
	if(ratioRange < 5) h2->SetAxisRange(0,5,"Y");
	else h2->SetAxisRange(0,1.2*ratioRange,"Y");
	h2->DrawCopy("le1");
	jumSun(xMin,1,xMax,1);
	c-> SaveAs(Form("pdf/compare_%s_cut_%s.pdf",var.Data(),theCut.GetTitle()));
}
