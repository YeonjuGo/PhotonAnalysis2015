// Author : Yeonju Go
// draw basic distributions (pthat, jtpt, refpt, photon pt, gen photon pt) simulataneously for the two different samples.
// two samples are forest(not skimmed) 
//
// leading photon condition X, only pthat30 
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TCut.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "stdio.h"
#include <iostream>

void compareDist1(){
//	TFile* myFile = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD30/JEC_pPbAllQCD30.root");
//	TFile* alexFile = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root");
//	TFile* myFile = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD50/JEC_pPbAllQCD50.root");
//	TFile* alexFile = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root");
//	TFile* myFile = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD80/JEC_pPbAllQCD80.root");
//      TFile* alexFile = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root");
	TFile* myFile = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD120/JEC_pPbAllQCD120.root");
        TFile* alexFile = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root");

	TTree* myTree = (TTree*) myFile -> Get("ak3PFJetAnalyzer/t");
	TTree* alexTree = (TTree*) alexFile -> Get("ak3PFJetAnalyzer/t");
	TTree* myTree_pho = (TTree*) myFile -> Get("multiPhotonAnalyzer/photon");
	TTree* alexTree_pho = (TTree*) alexFile -> Get("multiPhotonAnalyzer/photon");
	
	int myEntries = myTree->GetEntries();
	int alexEntries = alexTree->GetEntries();
	cout << "file #1 entries : " << myEntries << endl; 
	cout << "file #2 entries : " << alexEntries << endl; 
//=====================================
//pthat
//======================================

	TCanvas* c4 = new TCanvas();
	TH1D* mypthat = new TH1D("mypthat", "pthat;pthat(GeV);", 500,0,500);
	TH1D* alexpthat = (TH1D*)mypthat -> Clone("alexpthat");
		
	myTree -> Draw("pthat>>+mypthat");
	mypthat = (TH1D*)gDirectory->Get("mypthat");
	alexTree -> Draw("pthat>>+alexpthat");
	alexpthat = (TH1D*)gDirectory->Get("alexpthat");

//	mypthat -> Scale(1./mypthat->Integral("width"));
//	alexpthat -> Scale(1./alexpthat->Integral("width"));
	mypthat -> Scale(1./myEntries);
	alexpthat -> Scale(1./alexEntries);

	mypthat -> SetLineColor(1);
	mypthat -> SetMarkerColor(1);
	alexpthat -> SetLineColor(2);
	alexpthat -> SetMarkerColor(2);

	cout << " mypthat area : " << mypthat->Integral("width")<<endl;
	cout << " alexpthat area : " << alexpthat->Integral("width")<<endl;
	alexpthat -> Draw();
	mypthat -> Draw("same hist");

//=====================================
//photon pt
//======================================

	TCanvas* c1 = new TCanvas();
	TH1D* mypt = new TH1D("mypt", "photon pt;photon pt(GeV);", 500,0,500);
	TH1D* alexpt = (TH1D*)mypt -> Clone("alexpt");
		
	myTree_pho -> Draw("pt>>+mypt", "(pt<1000)");
	mypt = (TH1D*)gDirectory->Get("mypt");
	alexTree_pho -> Draw("pt>>+alexpt", "(pt<1000)");
	alexpt = (TH1D*)gDirectory->Get("alexpt");

//	mypt -> Scale(1./mypt->Integral("width"));
//	alexpt -> Scale(1./alexpt->Integral("width"));
	mypt -> Scale(1./myEntries);
	alexpt -> Scale(1./alexEntries);

	mypt -> SetLineColor(1);
	mypt -> SetMarkerColor(1);
	alexpt -> SetLineColor(2);
	alexpt -> SetMarkerColor(2);

	cout << " mypt area : " << mypt->Integral("width")<<endl;
	cout << " alexpt area : " << alexpt->Integral("width")<<endl;
	alexpt -> Draw();
	mypt -> Draw("same hist");

//=====================================
//jtpt
//======================================

	TCanvas* c2 = new TCanvas();
	TH1D* myjtpt = new TH1D("myjtpt", "jtpt;jtpt(GeV);", 500,0,500);
	TH1D* alexjtpt = (TH1D*)myjtpt -> Clone("alexjtpt");
		
	myTree -> Draw("jtpt>>+myjtpt", "(jtpt<1000)");
	myjtpt = (TH1D*)gDirectory->Get("myjtpt");
	alexTree -> Draw("jtpt>>+alexjtpt", "(jtpt<1000)");
	alexjtpt = (TH1D*)gDirectory->Get("alexjtpt");

//	alexjtpt -> Scale(myjtpt->Integral(myjtpt->FindBin(25),myjtpt->FindBin(500),"width")/alexjtpt->Integral("width"));
//	myjtpt -> Scale(1./mypt->Integral("width"));

//	myjtpt -> Scale(1./myjtpt->Integral("width"));
//	alexjtpt -> Scale(1./alexjtpt->Integral("width"));

	myjtpt -> Scale(1./myEntries);
	alexjtpt -> Scale(1./alexEntries);

	myjtpt -> SetLineColor(1);
	myjtpt -> SetMarkerColor(1);
	alexjtpt -> SetLineColor(2);
	alexjtpt -> SetMarkerColor(2);

	cout << " myjtpt area : " << myjtpt->Integral("width")<<endl;
	cout << " alexjtpt area : " << alexjtpt->Integral("width")<<endl;
	alexjtpt -> Draw();
	myjtpt -> Draw("same hist");


//=====================================
//refpt
//======================================

	TCanvas* c3 = new TCanvas();
	TH1D* myrefpt = new TH1D("myrefpt", "refpt;refpt(GeV);", 500,0,500);
	TH1D* alexrefpt = (TH1D*)myrefpt -> Clone("alexrefpt");
		
	myTree -> Draw("refpt>>+myrefpt", "(refpt<1000 && refpt>0)");
	myrefpt = (TH1D*)gDirectory->Get("myrefpt");
	alexTree -> Draw("refpt>>+alexrefpt", "(refpt<1000 && refpt>0)");
	alexrefpt = (TH1D*)gDirectory->Get("alexrefpt");

	//alexrefpt -> Scale(myrefpt->Integral(myrefpt->FindBin(20),myrefpt->FindBin(500),"width")/alexrefpt->Integral("width"));
	//myrefpt -> Scale(1./mypt->Integral("width"));
	//alexrefpt -> Scale(1./alexpt->Integral("width"));

	myrefpt -> Scale(1./myEntries);
	alexrefpt -> Scale(1./alexEntries);

	myrefpt -> SetLineColor(1);
	myrefpt -> SetMarkerColor(1);
	alexrefpt -> SetLineColor(2);
	alexrefpt -> SetMarkerColor(2);

	cout << " myrefpt area : " << myrefpt->Integral("width")<<endl;
	cout << " alexrefpt area : " << alexrefpt->Integral("width")<<endl;
	alexrefpt -> Draw();
	myrefpt -> Draw("same hist");

}
