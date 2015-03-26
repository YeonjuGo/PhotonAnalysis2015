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
#include <TVectorT.h>
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

static const long MAXTREESIZE = 10000000000;

void ljet_subljet_vector(bool isMCcut=1, int WhichDphi=0){
	//'WhichDphi = 0' means that delta phi(gamma, jet) is over than 7*pi/8.
	//'WhichDphi = 1' means that delta phi(gamma, jet) is over than 2*pi/3.
	//'WhichDphi = 2' means that delta phi(gamma, jet) is over than 1*pi/2.
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.04);

	const double ptbins[] = {15,20,30,40,50,60,70,80,90,100,110,120,130,140,160,180,200,240,280};
	const int nptbins = sizeof(ptbins)/sizeof(double) - 1;

	TFile* fmc;
	TTree* tgj_mc;
	TFile* fdata;
	TTree* tgj_data;
	if(WhichDphi==0) {
		fmc = new TFile("/u/user/goyeonju/files/yskimfiles/pA/merged_yskim_HiForest_pPb_MIX_AllQCDPhoton_akPu3PF_AfterResCorr_final_subjet.root");
		tgj_mc = (TTree*) fmc->Get("tgj");	
		fdata = new TFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_pPb_DATA_photon30trig_localJEC_v1_subjet.root");
		tgj_data = (TTree*) fdata->Get("tgj");
	} else if(WhichDphi==1) {
		fmc = new TFile("/u/user/goyeonju/files/yskimfiles/pA/merged_yskim_HiForest_pPb_MIX_AllQCDPhoton_akPu3PF_AfterResCorr_final_subjet_dphi23.root");
		tgj_mc = (TTree*) fmc->Get("tgj");	
		fdata = new TFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_pPb_DATA_photon30trig_localJEC_v1_subjet_dphi23.root");
		tgj_data = (TTree*) fdata->Get("tgj");	
	} else {
		fmc = new TFile("/u/user/goyeonju/files/yskimfiles/pA/merged_yskim_HiForest_pPb_MIX_AllQCDPhoton_akPu3PF_AfterResCorr_final_subjet_dphi12.root");
		tgj_mc = (TTree*) fmc->Get("tgj");	
		fdata = new TFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_pPb_DATA_photon30trig_localJEC_v1_subjet_dphi12.root");
		tgj_data = (TTree*) fdata->Get("tgj");	
	}
	
	const int nPhoCut =3;
	TCut phoEtaCut = "abs(photonEta)<1.44";
	TCut phoEtCut[nPhoCut] = {"photonEt>40","photonEt>60","photonEt>80"} ;
	TCut candidateCut = "sigmaIetaIeta<0.01";
	TCut isoCut = "ecalIso<4.2&&hcalIso<2.2&&trackIso<2.0";

	TCut mcCut;
	if(isMCcut==1) mcCut = "genIso<5&&abs(genMomId)<=22&&genPhotonEt>0"; 
	else mcCut = ""; 

	TCut njet1 = "njet==1";
	TCut njet2 = "njet==2";
	TCut totCut[2];
	TCut totCut_mc[2];
	totCut[0] = njet1 && phoEtCut[2] && phoEtaCut && candidateCut && isoCut;
	totCut[1] = njet2 && phoEtCut[2] && phoEtaCut && candidateCut && isoCut;
	totCut_mc[0] = njet1 && phoEtCut[2] && phoEtaCut && candidateCut && isoCut && mcCut;
	totCut_mc[1] = njet2 && phoEtCut[2] && phoEtaCut && candidateCut && isoCut && mcCut;
	TCut totNonjetCut = phoEtCut[2] && phoEtaCut && candidateCut && isoCut;

	const int nCdt = 3;
	TString label[nCdt] = { "njet=1", "njet=2", "njet=1 or 2"};
	TH1D *hpho_data[nCdt];
	TH1D *hpho_mc[nCdt];
	TH1D *hjpt_data[nCdt];
	TH1D *hjpt_mc[nCdt];
	TH1D *hratio_data[nCdt];
	TH1D *hratio_mc[nCdt];
	for(int i=0;i<nCdt;i++){
		hpho_data[i]=new TH1D(Form("hpho_data_%d",i),Form("%s;p_{T}^{#gamma};Normalized Entries",label[i].Data()), nptbins, ptbins);
		hpho_mc[i]=(TH1D*)hpho_data[i]->Clone(Form("hpho_mc_%d",i));
		hjpt_data[i]=new TH1D(Form("hjpt_data_%d",i),Form("%s;#Sigma(p_{T}^{jet});Normalized Entries",label[i].Data()), nptbins, ptbins);
		hjpt_mc[i]=(TH1D*)hjpt_data[i]->Clone(Form("hjpt_mc_%d",i));
		hratio_data[i]=new TH1D(Form("hratio_data_%d",i),Form("%s;#Sigma(p_{T}^{jet})/p_{T}^{#gamma};Normalized Entries",label[i].Data()),50,0,2);
		hratio_mc[i]=(TH1D*)hratio_data[i]->Clone(Form("hratio_mc_%d",i));
	}

//********************* photon pt dist ***********************
#if 0
	TCanvas* c_pho = new TCanvas("c_pho", "photon p_{T} distribution", 300, 300);
	TH1D *hpho_data_total;
	TH1D *hpho_mc_total;
	hpho_data_total=(TH1D*)hpho_data[0]->Clone("hpho_data_total");
	hpho_mc_total=(TH1D*)hpho_mc[0]->Clone("hpho_mc_total");
	tgj_data->Project(hpho_data_total->GetName(),"photonEt",totNonjetCut && "njet>0 && njet<3");
	tgj_mc->Project(hpho_mc_total->GetName(),"photonEt","ptHatWeight"*(totNonjetCut && mcCut && "njet>0 && njet<3"));
	hpho_data_total->SetMarkerStyle(20);
	hpho_mc_total->SetLineColor(4);
	hpho_data_total->Draw();		
	hpho_mc_total->Draw("same hist");		
#endif
	TString jetSum = "sqrt(lJetPt*lJetPt + sublJetPt*sublJetPt + 2*lJetPt*sublJetPt*cos(lJetPhi-sublJetPhi))";
	TString jetSumOverPho = "sqrt(lJetPt*lJetPt + sublJetPt*sublJetPt + 2*lJetPt*sublJetPt*cos(lJetPhi-sublJetPhi))/photonEt";
#if 0	
//********************* leading jet pt vs. sum of leading and sub-leading jet pt 2D ***********************
	//TString jetSum = "sqrt(lJetPt*lJetPt + sublJetPt*sublJetPt)";
	TCanvas* clsubl = new TCanvas("clsubl", "leading jet vs. sum jet pt (2D)", 500, 500);
	TCanvas* clsubl_mc = new TCanvas("clsubl_mc", "leading jet vs. sum jet pt (2D)", 500, 500);
	TH2D *hljet_sumjet_data=new TH2D("hljet_sumjet_data","hljet_sumjet data;leading jet pt;sum jet pt", nptbins, ptbins, nptbins, ptbins);
	TH2D *hljet_sumjet_mc=new TH2D("hljet_sumjet_mc","hljet_sumjet mc;leading jet pt;sum jet pt", nptbins, ptbins, nptbins, ptbins);
	tgj_data->Draw(Form("%s:lJetPt >> hljet_sumjet_data",jetSum.Data()),totNonjetCut && njet2);
	tgj_mc->Draw(Form("%s:lJetPt >> hljet_sumjet_mc",jetSum.Data()),"ptHatWeight"*(totNonjetCut && mcCut && njet2));
	clsubl->cd();
	hljet_sumjet_data->Draw("colz");		
	clsubl_mc->cd();
	hljet_sumjet_mc->Draw("colz");	
	clsubl->SaveAs("pdf/leading_sumjet_2D_data.pdf");
	clsubl_mc->SaveAs("pdf/leading_sumjet_2D_mc.pdf");
//********************* leading jet pt vs. sum of leading and sub-leading jet pt 2D ***********************
	TCanvas* c_2D_sub = new TCanvas("c2D_sub", "leading jet vs. sub-leading jet pt (2D)", 500, 500);
	TCanvas* c_2D_sub_mc = new TCanvas("c2D_sub_mc", "leading jet vs. sub-leading jet pt (2D)", 500, 500);
	TH2D *hljet_subljet_data=new TH2D("hljet_subljet_data","hljet_subljet data;leading jet pt;sub-leading jet pt", nptbins, ptbins, nptbins, ptbins);
	TH2D *hljet_subljet_mc=new TH2D("hljet_subljet_mc","hljet_subljet mc;leading jet pt;sub-leading jet pt", nptbins, ptbins, nptbins, ptbins);
	tgj_data->Draw("sublJetPt:lJetPt >> hljet_subljet_data",totNonjetCut && njet2);
	tgj_mc->Draw("sublJetPt:lJetPt >> hljet_subljet_mc","ptHatWeight"*(totNonjetCut && mcCut && njet2));
	c_2D_sub->cd();
	hljet_subljet_data->Draw("colz");		
	c_2D_sub_mc->cd();
	hljet_subljet_mc->Draw("colz");		
	c_2D_sub->SaveAs("pdf/leading_subleading_2D_data.pdf");
	c_2D_sub_mc->SaveAs("pdf/leading_subleading_2D_mc.pdf");
#endif

//========================================================================================================
//********************* hist filling  ***********************
//========================================================================================================
	TCanvas* ccc;
	tgj_data->Project(hpho_data[0]->GetName(),"photonEt",totCut[0]);
	tgj_mc->Project(hpho_mc[0]->GetName(),"photonEt","ptHatWeight"*(totCut_mc[0]));		
	tgj_data->Project(hpho_data[1]->GetName(),"photonEt",totCut[1]);
	tgj_mc->Project(hpho_mc[1]->GetName(),"photonEt","ptHatWeight"*(totCut_mc[1]));	
	tgj_data->Project(hjpt_data[0]->GetName(),"lJetPt",totCut[0]);
	tgj_mc->Project(hjpt_mc[0]->GetName(),"lJetPt","ptHatWeight"*(totCut_mc[0]));
	tgj_data->Project(hjpt_data[1]->GetName(),jetSum.Data(),totCut[1]);
	tgj_mc->Project(hjpt_mc[1]->GetName(),jetSum.Data(),"ptHatWeight"*(totCut_mc[1]));
	tgj_data->Project(hratio_data[0]->GetName(),jetSumOverPho.Data(),totCut[0]);
	tgj_mc->Project(hratio_mc[0]->GetName(),jetSumOverPho.Data(),"ptHatWeight"*(totCut_mc[0]));
	tgj_data->Project(hratio_data[1]->GetName(),jetSumOverPho.Data(),totCut[1]);
	tgj_mc->Project(hratio_mc[1]->GetName(),jetSumOverPho.Data(),"ptHatWeight"*(totCut_mc[1]));

	TH1D *temp_hpho_data[nCdt-1];
	TH1D *temp_hpho_mc[nCdt-1];
	TH1D *temp_hjpt_data[nCdt-1];
	TH1D *temp_hjpt_mc[nCdt-1];
	TH1D *temp_hratio_data[nCdt-1];
	TH1D *temp_hratio_mc[nCdt-1];
	
	for(int i=0;i<2;i++){
		temp_hpho_data[i]=(TH1D*)hjpt_data[i]->Clone(Form("temp_hjpt_data_%d",i));
		temp_hpho_mc[i]=(TH1D*)hjpt_mc[i]->Clone(Form("temp_hjpt_mc_%d",i));
		temp_hjpt_data[i]=(TH1D*)hjpt_data[i]->Clone(Form("temp_hjpt_data_%d",i));
		temp_hjpt_mc[i]=(TH1D*)hjpt_mc[i]->Clone(Form("temp_hjpt_mc_%d",i));
		temp_hratio_data[i]=(TH1D*)hratio_data[i]->Clone(Form("temp_hratio_data_%d",i));
		temp_hratio_mc[i]=(TH1D*)hratio_mc[i]->Clone(Form("temp_hratio_mc_%d",i));

		hpho_data[i]->Scale(1.,"width");
		hpho_mc[i]->Scale(1.,"width");
		hjpt_data[i]->Scale(1.,"width");
		hjpt_mc[i]->Scale(1.,"width");
		hratio_data[i]->Scale(1.,"width");
		hratio_mc[i]->Scale(1.,"width");

		temp_hpho_data[i]->SetMarkerStyle(20);
		temp_hjpt_data[i]->SetMarkerStyle(20);
		temp_hratio_data[i]->SetMarkerStyle(20);
		hpho_data[i]->SetMarkerStyle(20);
		hjpt_data[i]->SetMarkerStyle(20);
		hratio_data[i]->SetMarkerStyle(20);

		temp_hpho_mc[i]->SetLineColor(4);
		temp_hjpt_mc[i]->SetLineColor(4);
		temp_hratio_mc[i]->SetLineColor(4);
		hpho_mc[i]->SetLineColor(4);
		hjpt_mc[i]->SetLineColor(4);
		hratio_mc[i]->SetLineColor(4);
	}
	hpho_data[2]->Add(hpho_data[0],hpho_data[1]);
	hpho_mc[2]->Add(hpho_mc[0],hpho_mc[1]);
	hjpt_data[2]->Add(hjpt_data[0],hjpt_data[1]);
	hjpt_mc[2]->Add(hjpt_mc[0],hjpt_mc[1]);
	hratio_data[2]->Add(hratio_data[0],hratio_data[1]);
	hratio_mc[2]->Add(hratio_mc[0],hratio_mc[1]);
/*
	hpho_data[2]->Scale(1.,"width");
	hpho_mc[2]->Scale(1.,"width");
	hjpt_data[2]->Scale(1.,"width");
	hjpt_mc[2]->Scale(1.,"width");
	hratio_data[2]->Scale(1.,"width");
	hratio_mc[2]->Scale(1.,"width");
*/	
	hpho_data[2]->SetMarkerStyle(20);	
	hjpt_data[2]->SetMarkerStyle(20);
	hratio_data[2]->SetMarkerStyle(20);
	hpho_mc[2]->SetLineColor(4);
	hjpt_mc[2]->SetLineColor(4);
	hratio_mc[2]->SetLineColor(4);

	TLegend *l1 = new TLegend(0.5365615,0.6445304,0.8577623,0.846736,NULL,"brNDC");
	easyLeg(l1);
	l1->AddEntry(hjpt_data[0],"DATA");	
	l1->AddEntry(hjpt_mc[0],"MC");	
	TCanvas *c1 = new TCanvas("c1", "vector sum of jet pt dist", 1500,900);
	c1 -> Divide(3,3);

	c1->cd(1);
	hjpt_data[0]->Draw();
	hjpt_mc[0]->Draw("same hist");
	c1->cd(2);
	hjpt_data[1]->Draw();
	hjpt_mc[1]->Draw("same hist");
	c1->cd(3);
	hjpt_data[2]->Draw();
	hjpt_mc[2]->Draw("same hist");
	c1->cd(4);
	hratio_data[0]->Draw();
	hratio_mc[0]->Draw("same hist");
	c1->cd(5);
	hratio_data[1]->Draw();
	hratio_mc[1]->Draw("same hist");
	c1->cd(6);
	hratio_data[2]->Draw();
	hratio_mc[2]->Draw("same hist");

	c1->cd(7);
	hpho_data[0]->Draw();
	hpho_mc[0]->Draw("same hist");
	c1->cd(8);
	hpho_data[1]->Draw();
	hpho_mc[1]->Draw("same hist");
	c1->cd(9);
	hpho_data[2]->Draw();
	hpho_mc[2]->Draw("same hist");

	c1->SaveAs(Form("pdf/vectorSum_ratio_total_isMCcut%d_dphi%d.pdf",(int)isMCcut,WhichDphi));
#if 0
	TCanvas* can=new TCanvas("can", "", 1400,1350);
	can ->Divide(4,3);
	
	can->cd(1);
	temp_hjpt_data[0]->Draw();
	temp_hjpt_mc[0]->Draw("same hist");
	can->cd(2);
	temp_hjpt_data[1]->Draw();
	temp_hjpt_mc[1]->Draw("same hist");
	can->cd(3);
	temp_hratio_data[0]->Draw();
	temp_hratio_mc[0]->Draw("same hist");
	can->cd(4);
	temp_hratio_data[1]->Draw();
	temp_hratio_mc[1]->Draw("same hist");

	can->cd(5);
	hjpt_data[0]->Draw();
	hjpt_mc[0]->Draw("same hist");
	can->cd(6);
	hjpt_data[1]->Draw();
	hjpt_mc[1]->Draw("same hist");
	can->cd(7);
	hratio_data[0]->Draw();
	hratio_mc[0]->Draw("same hist");
	can->cd(8);
	hratio_data[1]->Draw();
	hratio_mc[1]->Draw("same hist");

	can->cd(9);
	hjpt_data[2]->Draw();
	hjpt_mc[2]->Draw("same hist");
	can->cd(10);
	hratio_data[2]->Draw();
	hratio_mc[2]->Draw("same hist");
#endif	
}

