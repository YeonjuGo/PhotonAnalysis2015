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

void spectrum_closure(bool doWeight=1){

	//TFile* myFile = new TFile("/u/user/goyeonju/2014/141015_JEC/CMSSW_5_3_20/src/JetMETAnalysis/JetAnalyzers/bin/merged_pPb_leadingpho_photonThreshold30.root");
	//TFile* alexFile = new TFile("/u/user/goyeonju/2014/141215_JECClosure/merged_pPb_leadingpho_eventselection_photonThreshold30.root");
	
	//only pthat 30 GeV fiils.
	//TFile* myFile = new TFile("/u/user/goyeonju/2014/141015_JEC/CMSSW_5_3_20/src/JetMETAnalysis/JetAnalyzers/bin/test/leadingpho40_pPb30.root");
	TFile* myFile = new TFile("test/pPb30_photonThreshold35_temp.root");
//	TFile* alexFile = new TFile("/u/user/goyeonju/2014/141215_JECClosure/test/leadingpho40_pPb30_photonThreshold25_temp.root");
	TFile* alexFile = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root");
	//TFile* alexFile = new TFile("/u/user/goyeonju/2014/141215_JECClosure/test/leadingpho40_pPb30_temp.root");
	TTree* myTree = (TTree*) myFile -> Get("ak3PFJetAnalyzer/t");
	TTree* alexTree = (TTree*) alexFile -> Get("ak3PFJetAnalyzer/t");

//=====================================
//photon pt
//======================================
/*
	TCanvas* c1 = new TCanvas();
	TH1D* mypt = new TH1D("mypt", "pt", 500,0,500);
	TH1D* alexpt = (TH1D*)mypt -> Clone("alexpt");
		
	if(doWeight==1)myTree -> Draw("pt>>+mypt", "weight*(pt<1000)");
	else myTree -> Draw("pt>>+mypt", "(pt<1000)");
	mypt = (TH1D*)gDirectory->Get("mypt");
	if(doWeight==1) alexTree -> Draw("pt>>+alexpt", "weight*(pt<1000)");
	else alexTree -> Draw("pt>>+alexpt", "(pt<1000)");
	alexpt = (TH1D*)gDirectory->Get("alexpt");

	mypt -> Scale(1./mypt->Integral("width"));
	alexpt -> Scale(1./alexpt->Integral("width"));

	mypt -> SetLineColor(1);
	mypt -> SetMarkerColor(1);
	alexpt -> SetLineColor(2);
	alexpt -> SetMarkerColor(2);

	cout << " mypt area : " << mypt->Integral("width")<<endl;
	cout << " alexpt area : " << alexpt->Integral("width")<<endl;
	alexpt -> Draw();
	mypt -> Draw("same hist");
*/
//=====================================
//jtpt
//======================================

	TCanvas* c2 = new TCanvas();
	TH1D* myjtpt = new TH1D("myjtpt", "jtpt", 500,0,500);
	TH1D* alexjtpt = (TH1D*)myjtpt -> Clone("alexjtpt");
		
	if(doWeight==1) myTree -> Draw("jtpt>>+myjtpt", "weight*(jtpt<1000)");
	else myTree -> Draw("jtpt>>+myjtpt", "(jtpt<1000)");
	myjtpt = (TH1D*)gDirectory->Get("myjtpt");
	if(doWeight==1) alexTree -> Draw("jtpt>>+alexjtpt", "weight*(jtpt<1000)");
	else alexTree -> Draw("jtpt>>+alexjtpt", "(jtpt<1000)");
	alexjtpt = (TH1D*)gDirectory->Get("alexjtpt");

	myjtpt -> Scale(1./myjtpt->Integral("width"));
	alexjtpt -> Scale(1./alexjtpt->Integral("width"));

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
	TH1D* myrefpt = new TH1D("myrefpt", "refpt", 500,0,500);
	TH1D* alexrefpt = (TH1D*)myrefpt -> Clone("alexrefpt");
		
	if(doWeight==1) myTree -> Draw("refpt>>+myrefpt", "weight*(refpt<1000 && refpt>0)");
	else myTree -> Draw("refpt>>+myrefpt", "(refpt<1000 && refpt>0)");
	myrefpt = (TH1D*)gDirectory->Get("myrefpt");
	if(doWeight==1) alexTree -> Draw("refpt>>+alexrefpt", "weight*(refpt<1000 && refpt>0)");
	else alexTree -> Draw("refpt>>+alexrefpt", "(refpt<1000 && refpt>0)");
	alexrefpt = (TH1D*)gDirectory->Get("alexrefpt");

	myrefpt -> Scale(1./myrefpt->Integral("width"));
	alexrefpt -> Scale(1./alexrefpt->Integral("width"));

	myrefpt -> SetLineColor(1);
	myrefpt -> SetMarkerColor(1);
	alexrefpt -> SetLineColor(2);
	alexrefpt -> SetMarkerColor(2);

	cout << " myrefpt area : " << myrefpt->Integral("width")<<endl;
	cout << " alexrefpt area : " << alexrefpt->Integral("width")<<endl;
	alexrefpt -> Draw();
	myrefpt -> Draw("same hist");

//=====================================
//pthat
//======================================

	TCanvas* c4 = new TCanvas();
	TH1D* mypthat = new TH1D("mypthat", "pthat", 500,0,500);
	TH1D* alexpthat = (TH1D*)mypthat -> Clone("alexpthat");
		
	if(doWeight==1) myTree -> Draw("pthat>>+mypthat", "weight*(pthat<1000)");
	else myTree -> Draw("pthat>>+mypthat", "(pthat<1000)");
	mypthat = (TH1D*)gDirectory->Get("mypthat");
	if(doWeight==1) alexTree -> Draw("pthat>>+alexpthat", "weight*(pthat<1000)");
	else alexTree -> Draw("pthat>>+alexpthat", "(pthat<1000)");
	alexpthat = (TH1D*)gDirectory->Get("alexpthat");

	mypthat -> Scale(1./mypthat->Integral("width"));
	alexpthat -> Scale(1./alexpthat->Integral("width"));

	mypthat -> SetLineColor(1);
	mypthat -> SetMarkerColor(1);
	alexpthat -> SetLineColor(2);
	alexpthat -> SetMarkerColor(2);

	cout << " mypthat area : " << mypthat->Integral("width")<<endl;
	cout << " alexpthat area : " << alexpthat->Integral("width")<<endl;
	alexpthat -> Draw();
	mypthat -> Draw("same hist");

}
