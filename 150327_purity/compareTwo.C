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

void compareTwo(){
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.05);

	const double ptbins[] = {15,20,30,40,50,60,70,80,90,100,110,120,130,140,160,180,200,240,280};
	const int nptbins = sizeof(ptbins)/sizeof(double) - 1;
	double AvePtBin[nptbins];

	for(int i=0;i<nptbins;i++){
		AvePtBin[i] = (ptbins[i+1]+ptbins[i])/2.0;
	}

	//old
	TFile * f1 = new TFile("/u/user/goyeonju/files/forest/pA/pA_photonSkimForest_v85_HLT_PAPhoton30_NoCaloIdVL_v1_highPtPhoton40.root");
	TTree* t1 = (TTree*) f1->Get("multiPhotonAnalyzer/photon");
	//new
	TFile * f2 = new TFile("/u/user/goyeonju/files/forest/pA/pPb_DATA_photon30trig_localJEC_v1.root");
	TTree* t2 = (TTree*) f2->Get("multiPhotonAnalyzer/photon");

	TH1D* h1 = new TH1D("h1", "", 100, 0, 1000);
	TH1D* h2 = new TH1D("h2", "", 100, 0, 1000);

	//configuration
	TString var = "pt";
	TCut cut = "isEle==0";

	//draw
	t1->Draw(Form("%s>>h1",var.Data()),cut );
	t2->Draw(Form("%s>>h2",var.Data()),cut );

	TCanvas* c1;
	h1->SetMarkerStyle(20);
	h1->Draw("p");
	h2->Draw("same hist");


//	c1 ->SaveAs(Form("pdf/compareTwo_var_%s_cut_%s.pdf",var.Data(),cut.GetTitle()));
}

