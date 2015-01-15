//Author : Yeonju Go

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
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
#include <sstream>

void FillMeanSigma(Int_t ip, TH1F *h1F, TH1F *hArM, TH1F *hRMS, TH1F *hMean, TH1F *hSigma);
TF1* FillGaussMeanSigma(Int_t ip, TH1F *h1F, TH1F *hMean, TH1F *hSigma);
//void FitGauss(TH1F* inHist_p, Float_t &mean, Float_t &meanError, Float_t &res, Float_t &resError);

const Int_t maxEntry = 5; //if fewer than this number of entries, ignore histogram
//const double fitmin=0.90;
//const double fitmax=1.10;
const double fitmin=0.00;
const double fitmax=2.0;
const TString fopt="RQ+";
const Int_t iFit=0;
const Int_t knpx=2000;

const Double_t ptbins[] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280};
//const Double_t ptbins[] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,340, 400, 500};
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
  
void leadingPho_draw_JEC(
        const char* infile="/home/goyeonju/CMS/2015/gammaJetAnalysis/histogramProducer/forestFiles/pA/PA2013_pyquen_allQCDPhoton30to50_forestv85.root", 
        TString calgo="ak3PF", bool savePlots=1, bool gausfitting=1)
{

    TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetMarkerStyle(20);

	TFile* fin = TFile::Open(infile,"readonly");
    cout << infile << endl;
    
	TTree* tree = (TTree*) fin -> Get(calgo+"JetAnalyzer/t");
	//TTree* tree = (TTree*) fin -> Get(calgo+"JetAnalyzer/t"Form("%sJetAnalyzer/t", calgo.str().c_str())); 
    TString calgoDir;
    if(gausfitting==1) calgoDir =calgo+"_gaus";
    if(gausfitting==0) calgoDir =calgo+"_mean";

//===========================================================================
// to test weighting factor
//===========================================================================
/*	TH1D *histpthat = new TH1D("histpthat",";pthat(GeV);arbitrary entries",500,0,1000);
	tree->Draw("pthat>>histpthat");
	//tree->Draw("pthat>>histpthat","weight");
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
	TF1* fratio[nptbins];
	for(int i=0; i<nptbins; i++){
		TCut ptCut = Form("jtpt > %lf && jtpt < %lf",ptbins[i], ptbins[i+1]);
		TString hName = Form("reco_over_gen_pt_%d",(Int_t)ptbins[i]);
		TH1F* ratio;
		if(i==nptbins-1) ratio = new TH1F(hName, Form("%d GeV< gen p_{T}<500 GeV;reco/gen p_{T}",(Int_t)ptbins[i]), 200, fitmin, fitmax);
		else ratio = new TH1F(hName, Form("%d GeV< gen p_{T}<%d GeV;reco/gen p_{T}",(Int_t)ptbins[i],(Int_t)ptbins[i+1]), 200, fitmin, fitmax);

		//////////////// without weighting
		//tree -> Draw(Form("jtpt/refpt >> reco_over_gen_pt_%d",(Int_t)ptbins[i]), Form("(jtpt > %lf && jtpt < %lf && refpt>0 && refpt<1000 && jteta>-3.0 && jteta<3.0)",ptbins[i], ptbins[i+1]));
		//////////////// x axis is jtpt(reco pt)
		//tree -> Draw(Form("jtpt/refpt >> reco_over_gen_pt_%d",(Int_t)ptbins[i]), Form("weight*(jtpt > %lf && jtpt < %lf && refpt>0 && refpt<1000 && jteta>-3.0 && jteta<3.0)",ptbins[i], ptbins[i+1]));
		//////////////// x axis is refpt(gen pt)
		tree -> Draw(Form("jtpt/refpt >> reco_over_gen_pt_%d",(Int_t)ptbins[i]), Form("weight*(refpt > %lf && refpt < %lf && jtpt<1000 && jteta>-3.0 && jteta<3.0)",ptbins[i], ptbins[i+1]));
		//////////////// HFsum > 20 
	//	tree -> Draw(Form("jtpt/refpt >> reco_over_gen_pt_%d",(Int_t)ptbins[i]), Form("weight*(jtpt > %lf && jtpt < %lf && refpt>0 && refpt<1000 && jteta>-3.0 && jteta<3.0 && hiHF>20)",ptbins[i], ptbins[i+1]));
		//////////////// HFsum < 20
		//tree -> Draw(Form("jtpt/refpt >> reco_over_gen_pt_%d",(Int_t)ptbins[i]), Form("weight*(jtpt > %lf && jtpt < %lf && refpt>0 && refpt<1000 && jteta>-3.0 && jteta<3.0 && hiHF<=20)",ptbins[i], ptbins[i+1]));
		
		ratio = (TH1F*)gDirectory->Get(hName);
		if(gausfitting==0)	FillMeanSigma(i, ratio, hArM, hRMS, hMean, hSigma);
		else if(gausfitting==1) fratio[i]= FillGaussMeanSigma(i, ratio, hMean, hSigma);
    #if 0
		d[i] = new TCanvas(Form("reco/gen_%d",(Int_t)ptbins[i]),Form("reco/gen_%d",(Int_t)ptbins[i]));
		ratio->DrawClone();
		if(savePlots)
			d[i]->SaveAs(calgoDir+"/"+hName+".gif");
			//d[i]->SaveAs(Form("%s/"+hName+".gif",calgoDir));
    #endif
	}

    TLegend*l1 = new TLegend(0.65, 0.60, 0.85, 0.80, calgoDir);
	TCanvas *c[4];
	c[0] = new TCanvas("c_hMean","c_hMean");
	c[0] -> cd();
	//hMean->SetAxisRange(0.99,1.10,"Y");
	hMean->Draw();
	l1->Draw();
	c[1] = new TCanvas("c_hSigma","c_hSigma");
	hSigma->SetAxisRange(0.0,0.30,"Y");
	hSigma->Draw();
	l1->Draw();
	if(savePlots)
	{
		c[0]->SaveAs(calgoDir+"/hMean.gif");
		c[1]->SaveAs(calgoDir+"/hSigma.gif");
	    //c[0]->SaveAs(Form("%s/hMean.gif",calgoDir));
		//c[1]->SaveAs(Form("%s/hSigma.gif",calgoDir));
	}    

}
#if 1
void FillMeanSigma(Int_t ip,TH1F *h1F,TH1F *hArM,TH1F *hRMS,TH1F *hMean,TH1F *hSigma)
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
#endif

TF1* FillGaussMeanSigma(Int_t ip, TH1F *h1F, TH1F *hMean, TH1F *hSigma)
{
	TF1* myFit = new TF1("myFit", "gaus");
	//myFit -> SetParameters(h1F->GetMaximum(),1.0, 0.1);
	h1F->Fit(myFit->GetName(), "QL", "");
#if 1	
	if(h1F->Integral()>0){
		hMean ->SetBinContent  (ip+1,h1F->GetFunction(myFit->GetName())->GetParameter(1));
		hMean ->SetBinError  (ip+1,h1F->GetFunction(myFit->GetName())->GetParError(1));
		hSigma->SetBinContent  (ip+1,h1F->GetFunction(myFit->GetName())->GetParameter(2));
		hSigma->SetBinError  (ip+1,h1F->GetFunction(myFit->GetName())->GetParError(2));
	}
#endif

	return myFit;
}

#if 0
void FillGaussMeanSigma(Int_t ip, TH1F *h1F2, TH1F *hMean2, TH1F *hSigma2)
{
//	TCanvas* c1 = new TCanvas();
//	c1 -> cd();
//	h1F -> Draw();
#if 0	h1F->Fit("gaus", "Q L", "");
	if(h1F->Integral()>0){
		hMean ->SetBinContent  (ip+1,h1F->GetFunction("gaus")->GetParameter(1));
		hMean ->SetBinError  (ip+1,h1F->GetFunction("gaus")->GetParError(1));
		hSigma->SetBinContent  (ip+1,h1F->GetFunction("gaus")->GetParameter(2));
		hSigma->SetBinError  (ip+1,h1F->GetFunction("gaus")->GetParError(2));
	}
#endif
}
#endif
#if 0
void FitGauss(TH1F* inHist_p, Float_t &mean, Float_t &meanError, Float_t &res, Float_t &resError)
{
  inHist_p->Fit("gaus", "Q L", "");

  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);

  Float_t prob = inHist_p->GetFunction("gaus")->GetProb();
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

  inHist_p->Fit("gaus", "Q L", "", minPt, maxPt);

  if(TMath::Abs(1.00 - mean) < TMath::Abs(1.00 - inHist_p->GetFunction("gaus")->GetParameter(1)) && prob > 0.0001){
    inHist_p->Fit("gaus", "Q L", "");
    return;
  }

  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);

  return;
}
#endif
