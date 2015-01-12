//Author : Yeonju Go

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

void FillMeanSigma(Int_t ip, TH1 *h1F, TH1 *hArM, TH1 *hRMS, TH1 *hMean, TH1 *hSigma);

const Int_t maxEntry = 5; //if fewer than this number of entries, ignore histogram
const double fitmin=0.00;
const double fitmax=3.00;
const TString fopt="RQ+";
const Int_t iFit=0;
const Int_t knpx=2000;

const Double_t ptbins[] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,340, 400, 500};
//double ptbins[] ={10,15,20,27,35,45,57,72,90,120,150,200,300,400,550,750,1000};
const int nptbins = sizeof(ptbins)/sizeof(double) - 1;

const int filepthat[]={30,50,80,120,170,9999};
const int nFile = 5;

//for 2013 HiWinter pPb
//const int kAlgos = 12; 
//const char *calgo[kAlgos] = {"ak3PF","ak4PF","ak5PF","ak3Calo","ak4Calo","ak5Calo",
//                     "akPu3PF","akPu4PF","akPu5PF","akPu3Calo","akPu4Calo","akPu5Calo"
//};

//const double weight[]={1., 1., 1., 1., 1.};
const double weight[]={102400.0/102400.0, 39656.0/104640.0, 10157.0/104640.0, 2517.0/104640.0, 649.0/104640.0};
  
void draw_JEC(const char *calgo="ak3PF", bool savePlots=1){
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetMarkerStyle(20);
	//gStyle->SetTitleXSize(28);

	TFile* fin[5];
	fin[0] = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton30.root");
	fin[1] = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton50.root");
	fin[2] = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton80.root");
	fin[3] = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton120.root");
	fin[4] = new TFile("/u/user/goyeonju/files/forest/pA/forJEC/HiForest_pPb_MIX_AllQCDPhoton170.root");

	TTree* tree[5];
	tree[0] = (TTree*) fin[0] -> Get(Form("%sJetAnalyzer/t", calgo)); 
	tree[1] = (TTree*) fin[1] -> Get(Form("%sJetAnalyzer/t", calgo)); 
	tree[2] = (TTree*) fin[2] -> Get(Form("%sJetAnalyzer/t", calgo)); 
	tree[3] = (TTree*) fin[3] -> Get(Form("%sJetAnalyzer/t", calgo)); 
	tree[4] = (TTree*) fin[4] -> Get(Form("%sJetAnalyzer/t", calgo)); 

//===========================================================================
// to test weighting factor
//===========================================================================
/*	TH1D *histpthat = new TH1D("histpthat",";pthat(GeV);arbitrary entries",500,0,500);
	for(int ipt=0;ipt<5;ipt++){
		TH1D *htmp = (TH1D*)histpthat->Clone();
		htmp->SetName(Form("my_htmp_%d",ipt));
		htmp->Reset();
		tree[ipt]->Draw(Form("%s>>my_htmp_%d","pthat",ipt),Form("pthat>=%d && pthat<%d",filepthat[ipt], filepthat[ipt+1]));
		htmp->Scale(weight[ipt]);
		histpthat->Add(htmp);
		delete htmp;
	} 
	TCanvas *cpthat = new TCanvas("cpthat", "cpthat", 600, 600);
	histpthat -> SetMarkerStyle(20);
	histpthat -> Draw();
*/
//===========================================================================
// to derive JEC
//===========================================================================
	TH1F *hMean = new TH1F("hMean",";gen p_{T}(GeV);<reco p_{T}/gen p_{T}>",nptbins,ptbins);
	TH1F *hArM = new TH1F("hArM",";gen p_{T}(GeV);hArM <reco p_{T}/gen p_{T}>",nptbins,ptbins);
	TH1F *hSigma = new TH1F("hSigma",";gen p_{T}(GeV);#sigma(reco p_{T}/gen p_{T})",nptbins,ptbins);  
	TH1F *hRMS = new TH1F("hRMS",";gen p_{T}(GeV);RMS(reco p_{T}/gen p_{T})",nptbins,ptbins);
	
	TCut analysisCut = "jteta>-3.0 && jteta<3.0";

	TCanvas *d[nptbins];
	for(int i=0; i<nptbins; i++){
		TCut ptCut = Form("refpt > %lf && refpt < %lf",ptbins[i], ptbins[i+1]);
		TString hName = Form("reco_over_gen_pt_%d",(Int_t)ptbins[i]);
		TH1D* ratio;
		if(i==nptbins-1) ratio = new TH1D(hName, Form("%d GeV< gen p_{T}<500 GeV;reco/gen p_{T}",(Int_t)ptbins[i]), 200, fitmin, fitmax);
		else ratio = new TH1D(hName, Form("%d GeV< gen p_{T}<%d GeV;reco/gen p_{T}",(Int_t)ptbins[i],(Int_t)ptbins[i+1]), 200, fitmin, fitmax);

		for(int ipt=0;ipt<5;ipt++){
			TH1D *htmp = (TH1D*)ratio->Clone();
			htmp->SetName(Form("my_htmp_%d",ipt));
			htmp->Reset();
			TCut pthatCut = Form("pthat>=%d && pthat<%d",filepthat[ipt], filepthat[ipt+1]);
			tree[ipt]->Draw(Form("%s>>my_htmp_%d","jtpt/refpt",ipt), analysisCut&&ptCut&&pthatCut);
			htmp->Scale(weight[ipt]);
			ratio->Add(htmp);
			delete htmp;
		}        

		FillMeanSigma(i, ratio, hArM, hRMS, hMean, hSigma);

/*		d[i] = new TCanvas(Form("reco/gen_%d",(Int_t)ptbins[i]),Form("reco/gen_%d",(Int_t)ptbins[i]));
		ratio->DrawClone();
		if(savePlots)
			d[i]->SaveAs(Form("%s/"+hName+".gif",calgo));
*/	}

	TLegend *l1 = new TLegend(0.65, 0.60, 0.85, 0.80, calgo);
	TCanvas *c[4];
	c[0] = new TCanvas("c_hMean","c_hMean");
	hMean->SetAxisRange(0.99,1.10,"Y");
	hMean->Draw();
	l1->Draw();
	c[1] = new TCanvas("c_hSigma","c_hSigma");
	hSigma->SetAxisRange(0.0,0.30,"Y");
	hSigma->Draw();
	l1->Draw();
	if(savePlots)
	{
		c[0]->SaveAs(Form("%s/hMean.gif",calgo));
		c[1]->SaveAs(Form("%s/hSigma.gif",calgo));
	}    
}

void FillMeanSigma(Int_t ip,TH1 *h1F,TH1 *hArM,TH1 *hRMS,TH1 *hMean,TH1 *hSigma)
{
	if(h1F->GetEntries()<maxEntry){
		h1F->Scale(0.);
		hArM  ->SetBinContent(ip+1,-9999);
		hArM  ->SetBinError  (ip+1,0);
		hRMS  ->SetBinContent(ip+1,-9999);
		hRMS  ->SetBinError  (ip+1,0);
	}
	if(h1F->Integral()>0){
		hMean ->SetBinContent  (ip+1,h1F->GetMean());
		hMean ->SetBinError  (ip+1,h1F->GetMeanError());
		hSigma->SetBinContent  (ip+1,h1F->GetRMS());
		hSigma->SetBinError  (ip+1,h1F->GetRMSError());
	}
}
