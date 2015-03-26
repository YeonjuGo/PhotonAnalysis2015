// Author : Yeonju Go
// 20150316 incomplete
//
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
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>
//random 
#include <TRandom.h>
#include <TStopwatch.h>
#include <ctime>

static const long MAXTREESIZE = 10000000000;

void ljet_subljet(){
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.04);

	const double ptbins[] = {15,20,30,40,50,60,70,80,90,100,110,120,130,140,160,180,200,240,280};
	const int nptbins = sizeof(ptbins)/sizeof(double) - 1;

	TFile* fmc = new TFile("/u/user/goyeonju/files/yskimfiles/pA/merged_yskim_HiForest_pPb_MIX_AllQCDPhoton_akPu3PF_AfterResCorr_final_subjet.root");
	TTree* tgj_mc = (TTree*) fmc->Get("tgj");

	TFile* fdata = new TFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_pPb_DATA_photon30trig_localJEC_v1_subjet.root");
	TTree* tgj_data = (TTree*) fdata->Get("tgj");

	const int nPhoCut =3;
	TCut phoEtaCut = "abs(photonEta)<1.44";
	TCut phoEtCut[nPhoCut] = {"photonEt>40","photonEt>60","photonEt>80"} ;
	TCut candidateCut = "sigmaIetaIeta<0.01";
	TCut isoCut = "ecalIso<4.2&&hcalIso<2.2&&trackIso<2.0";

	TCut mcCut = "genIso<5&&abs(genMomId)<=22&&genPhotonEt>0"; 
	TCut totCut = phoEtCut[2] && phoEtaCut && candidateCut && isoCut;

	TCanvas* c_temp = new TCanvas("c_temp", "jet p_{T} distribution", 300, 300);

	TH1D *hl_data = new TH1D("hl_data",";leading p_{T}^{Jet};Normalized Entries", nptbins,ptbins);
	TH1D *hl_mc = (TH1D*) hl_data->Clone("hl_mc");
	TH1D *hsubl_data = new TH1D("hsubl_data",";subleading p_{T}^{Jet};Normalized Entries", nptbins,ptbins);
	TH1D *hsubl_mc = (TH1D*) hsubl_data->Clone("hsubl_mc");

	tgj_data->Project(hl_data->GetName(),"lJetPt",totCut);
	tgj_mc->Project(hl_mc->GetName(),"lJetPt", Form("( %s ) * (%s && %s)","ptHatWeight", totCut.GetTitle(), mcCut.GetTitle()) );
	tgj_data->Project(hsubl_data->GetName(),"sublJetPt",totCut);
	tgj_mc->Project(hsubl_mc->GetName(),"sublJetPt", Form("( %s ) * (%s && %s)","ptHatWeight", totCut.GetTitle(), mcCut.GetTitle()) );

	hl_data->Scale(1./hl_data->Integral(),"width");
	hl_mc->Scale(1./hl_mc->Integral(),"width");
	hsubl_data->Scale(1./hsubl_data->Integral(),"width");
	hsubl_mc->Scale(1./hsubl_mc->Integral(),"width");

	TCanvas* c_cut = new TCanvas("c_cut", "jet p_{T} distribution", 800, 400); 
	c_cut->Divide(2,1);

	TLegend *l1 = new TLegend(0.5365615,0.6445304,0.8577623,0.846736,NULL,"brNDC");
	easyLeg(l1);
	//leading jet hist cosmetics
	c_cut->cd(1);
	double range = cleverRange(hl_data,hl_mc);
	hl_data->GetYaxis()->SetRangeUser(0.000001,range);
	hl_mc->GetYaxis()->SetRangeUser(0.000001,range);
	hl_data->SetMarkerStyle(20);
	hl_mc->SetLineColor(4);
	hl_mc->Draw("hist e");
	hl_data->Draw("same p");
	l1->AddEntry(hl_data,"Data","p");
	l1->AddEntry(hl_mc,"MC","l");
	l1->Draw("same");
	//drawText(str[j],0.6,0.75);			

	//subleading jet hist cosmetics
	c_cut->cd(2);
	range = cleverRange(hsubl_data,hsubl_mc);
	hsubl_data->GetYaxis()->SetRangeUser(0.000001,range);
	hsubl_mc->GetYaxis()->SetRangeUser(0.000001,range);
	hsubl_data->SetMarkerStyle(20);
	hsubl_mc->SetLineColor(4);
	hsubl_mc->Draw("hist e");
	hsubl_data->Draw("same p");
	c_cut ->SaveAs("pdf/leading_subleading_dist.pdf");

	//==============================================================
	// njet genMomId
	//==============================================================

	TH1D* hnjet_tot = new TH1D("hnjet_tot", "nJet distribution", 4,0,4);
	TH1D* hnjet_frag = (TH1D*) hnjet_tot->Clone("hnjet_frag");
	TH1D* hnjet_dir = (TH1D*) hnjet_tot->Clone("hnjet_dir");

	TCut dirCut = "abs(genMomId)==22";
	TCut fragCut = "abs(genMomId)<22";

	tgj_mc->Project(hnjet_tot->GetName(),"njet", Form("( %s ) * (%s && %s)","ptHatWeight", totCut.GetTitle(), mcCut.GetTitle()) );
	tgj_mc->Project(hnjet_frag->GetName(),"njet", Form("( %s ) * (%s && %s && %s)","ptHatWeight", totCut.GetTitle(), mcCut.GetTitle(), fragCut.GetTitle()) );
	tgj_mc->Project(hnjet_dir->GetName(),"njet", Form("( %s ) * (%s && %s&& %s)","ptHatWeight", totCut.GetTitle(), mcCut.GetTitle(), dirCut.GetTitle()) );

	TCanvas* c_njet = new TCanvas("c_njet", "# of jet distribution; # of jets;", 500, 500); 
	hnjet_dir->SetLineColor(2);
	hnjet_frag->SetLineColor(4);
	hnjet_tot->Draw("hist");
	hnjet_dir->Draw("same hist");
	hnjet_frag->Draw("same hist");
	gPad->SetLogy();	
	c_njet->SaveAs("pdf/njet_dist.pdf");

	//==============================================================
	// Evt size difference
	//==============================================================

	const Int_t CENTBINS[] = {0, 10, 20, 30, 200}; 
	const Int_t nCENTBINS = sizeof(CENTBINS)/sizeof(Int_t) -1;

	TH1D *hcentl_data[nCENTBINS];
	TH1D *hcentl_mc[nCENTBINS];
	TH1D *hsubcentl_data[nCENTBINS];
	TH1D *hsubcentl_mc[nCENTBINS];
	TCanvas* c_cent[nCENTBINS];

	for(Int_t j = 0; j < nCENTBINS; ++j) {
		c_cent[j] = new TCanvas(Form("c_cent%d",j), "jet p_{T} distribution", 800, 400); 
		c_cent[j]->Divide(2,1);

		hcentl_data[j] = new TH1D(Form("hcentl_data_cent%d",j),Form("%d<hf4Sum<%d;leading p_{T}^{Jet};Normalized Entries",CENTBINS[j],CENTBINS[j+1]), nptbins,ptbins);
		hcentl_mc[j]=(TH1D*) hcentl_data[j]->Clone(Form("hcentl_mc_cent%d",j));
		hsubcentl_data[j] = new TH1D(Form("hsubcentl_data_cent%d",j),Form("%d<hf4Sum<%d;subleading p_{T}^{Jet};Normalized Entries",CENTBINS[j],CENTBINS[j+1]), nptbins,ptbins);
		hsubcentl_mc[j]=(TH1D*) hsubcentl_data[j]->Clone(Form("hsubcentl_mc_cent%d",j));
	}

	for(Int_t j = 0; j < nCENTBINS; ++j) {
		TCut centCut = Form("((hf4Sum) >= %i) && ((hf4Sum) < %i)", CENTBINS[j], CENTBINS[j+1]);

		tgj_data->Project(hcentl_data[j]->GetName(),"lJetPt",totCut && centCut);
		tgj_mc->Project(hcentl_mc[j]->GetName(),"lJetPt", Form("( %s ) * (%s && %s && %s)","ptHatWeight", totCut.GetTitle(), mcCut.GetTitle(), centCut.GetTitle()) );
		tgj_data->Project(hsubcentl_data[j]->GetName(),"sublJetPt",totCut && centCut);
		tgj_mc->Project(hsubcentl_mc[j]->GetName(),"sublJetPt", Form("( %s ) * (%s && %s && %s)","ptHatWeight", totCut.GetTitle(), mcCut.GetTitle(), centCut.GetTitle()) );

		hcentl_data[j]->Scale(1./hcentl_data[j]->Integral(),"width");
		hcentl_mc[j]->Scale(1./hcentl_mc[j]->Integral(),"width");
		hsubcentl_data[j]->Scale(1./hsubcentl_data[j]->Integral(),"width");
		hsubcentl_mc[j]->Scale(1./hsubcentl_mc[j]->Integral(),"width");

		//hcentl_data[j]->GetYaxis()->SetOffset(0.000001,0.015);
		hcentl_data[j]->GetYaxis()->SetRangeUser(0.000001,0.015);
		hcentl_mc[j]->GetYaxis()->SetRangeUser(0.000001,0.015);
		hsubcentl_data[j]->GetYaxis()->SetRangeUser(0.000001,0.07);
		hsubcentl_mc[j]->GetYaxis()->SetRangeUser(0.000001,0.07);
		hcentl_data[j]->SetMarkerStyle(20);
		hcentl_mc[j]->SetLineColor(4);
		hsubcentl_data[j]->SetMarkerStyle(20);
		hsubcentl_mc[j]->SetLineColor(4);
	}
//	TLegend *l1 = new TLegend(0.5365615,0.6445304,0.8577623,0.846736,NULL,"brNDC");
///	easyLeg(l1);

	for(Int_t j = 0; j < nCENTBINS; ++j) {
		c_cent[j]->cd(1);

		hcentl_mc[j]->Draw("hist e");
		hcentl_data[j]->Draw("same p");
//		l1->AddEntry(hcentl_data[j],"Data","p");
//		l1->AddEntry(hcentl_mc[j],"MC","l");
		l1->Draw("same");

		//subleading jet hist cosmetics
		c_cent[j]->cd(2);
		hsubcentl_mc[j]->Draw("hist e");
		hsubcentl_data[j]->Draw("same p");
		c_cent[j]->SaveAs(Form("pdf/leading_subleading_dist_cent%d.pdf",j));

	}
}

