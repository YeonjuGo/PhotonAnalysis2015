// Author Yeonju Go
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TPad.h"
#include "stdio.h"
#include "gedPhotonUtility.h"
#include <iostream>
using namespace std;
//#include "../HiForestAnalysis/hiForest.h"

void evtSelEff()
{
	const TCut runCut = "run==181611";
	const TCut lumiCut = "lumi>=1 && lumi<=895";
	const TCut hltCut = "HLT_HIMinBiasHfOrBSC_v1==1";
	const TCut eventCut = hltCut;
	const int NFILE = 3; 
	const int NFILTER = 5;
	TCut evtFilter[NFILTER] = {"", "phfCoincFilter3==1","phfCoincFilter==1","pprimaryVertexFilter==1","pcollisionEventSelection==1"};
	TString evtFilterName[NFILTER] = {"NoFilter", "phfCoincFilter3","phfCoincFilter","pprimaryVertexFilter","pcollisionEventSelection"};
	int col[NFILTER]={kAzure,kGreen+1,kViolet};
	//int col[NFILTER]={kBlack,kAzure,kAzure-1,kGreen+1,kViolet};
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);

	TFile *f[NFILE];
	f[0]= new TFile("/u/user/goyeonju/files/centrality/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611_CMSSW5320_byYJ.root"); // data 53X
	f[1]= new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_7_5_0_pre5/src/centrality/crab_RECO/results/PbPb_minbias_data_2760_HIRun2011-25Aug2015_CMSSW750pre5_byYJ.root");
	f[2]= new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_7_5_0_pre5/src/centrality/crab_AOD/results/PbPb_minbias_data_2760_HIRun2011-25Aug2015_CMSSW750pre5_byYJ_AOD.root");

	TTree *t_evt[NFILE];
	TTree *t_skim[NFILE];
	TTree *t_hlt[NFILE];
	for(int i=0;i<NFILE;i++){
		t_evt[i] = (TTree*) f[i]->Get("hiEvtAnalyzer/HiTree");
		t_skim[i] = (TTree*) f[i]->Get("skimanalysis/HltTree");
		t_hlt[i] = (TTree*) f[i]->Get("hltanalysis/HltTree");
		t_evt[i]->AddFriend(t_skim[i]);
		t_evt[i]->AddFriend(t_hlt[i]);
	}

	TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
	//=========================
	// hiBin
	//=========================

	TH1D *h_hiBin[NFILE][NFILTER];
	for(int j=0; j<NFILE;j++){
		for(int i=0; i<NFILTER; i++)
		{
			h_hiBin[j][i] = new TH1D(Form("h_hiBin_ifile%d_filter%d",j,i), ";hiBin;",100,0,200);
			t_evt[j]->Draw(Form("hiBin>>+ %s",h_hiBin[j][i]->GetName()), eventCut && evtFilter[i]);
			h_hiBin[j][i] = (TH1D*)gDirectory->Get(h_hiBin[j][i]->GetName());
			h_hiBin[j][i] -> SetAxisRange(0,0.025,"Y");
			h_hiBin[j][i] -> GetYaxis()->SetRange(0,0.025);
			h_hiBin[j][i] -> GetYaxis()->SetRangeUser(0,0.025);
		}
	}

	cout << "End of fill histograms" << endl;
	TLegend* l1 = new TLegend(0.4, 0.65, 0.85, 0.80);
	l1 -> AddEntry(h_hiBin[0][0], "5_3_20 RECO");
	l1 -> AddEntry(h_hiBin[1][0], "7_5_0_pre5 RECO");
	l1 -> AddEntry(h_hiBin[2][0], "7_5_0_pre5 AOD");


	TH1D *h_temp = (TH1D*) h_hiBin[0][0]->Clone("h_temp");
	h_temp->SetAxisRange(0,0.025,"Y");
	TCanvas *c_hiBin[NFILTER];
	float xp=0.57;
	float yp=0.9;
	float dy=0.06;
	for(int i=0; i<NFILTER; i++){
		c_hiBin[i] = new TCanvas(Form("c_hiBin_filter%d",i), "hiBin", 400,800);
		c_hiBin[i]->Divide(1,2);
		c_hiBin[i]->cd(1);
		h_temp->DrawCopy();
		for(int j=0; j<NFILE;j++)
		{
			h_hiBin[j][i] -> SetMarkerStyle(20); 
			h_hiBin[j][i] -> SetMarkerSize(0.8); 
			h_hiBin[j][i] -> SetMarkerColor(col[j]); 
			h_hiBin[j][i] -> Scale(1./t_evt[j]->GetEntries(hltCut));
			if(j==0) h_hiBin[j][i] -> Draw("e1p");
			else h_hiBin[j][i] -> Draw("same && e1p");
		}
		l1->Draw("same");
		drawText(evtFilterName[i],xp,yp-0*dy);
	}

	TEfficiency* eff_hiBin[NFILE][NFILTER];
	TH1D* h_eff_hiBin[NFILE][NFILTER];
	for(int i=1; i<NFILTER; i++){
		for(int j=0; j<NFILE;j++){
			c_temp->cd();
			h_eff_hiBin[j][i]=(TH1D*) h_hiBin[j][0]->Clone(Form("h_eff_hiBin_ifile%d_ifilter%d",j,i));
			h_eff_hiBin[j][i]->SetTitle(";hiBin;Efficiency");
			h_eff_hiBin[j][i]->Divide(h_hiBin[j][i],h_hiBin[j][0],1,1,"b");
			h_eff_hiBin[j][i] -> SetMarkerStyle(20);
			h_eff_hiBin[j][i] -> SetMarkerSize(0.8);
			h_eff_hiBin[j][i] -> SetMarkerColor(col[j]);	
			h_eff_hiBin[j][i] -> SetAxisRange(0,1.0,"Y");
			c_hiBin[i]->cd(2);
			if(j==0) h_eff_hiBin[j][i]->Draw("e1p");
			else h_eff_hiBin[j][i]->Draw("same && e1p");
		}
		c_hiBin[i]->SaveAs(Form("pdf/hiBin_%s.png",evtFilterName[i].Data()));
	}

	//======================================
	//HF sum!!
	//======================================
#if 0

	TLine* t1 = new TLine(0,1,1000,1);
	t1->SetLineWidth(1);
	t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
	t1->SetLineColor(1); // 2 is red

	TH1D *HFsum_data[5];
	for(int i=0; i<Ncut; i++)
	{
		HFsum_data[i] = new TH1D(Form("HFsum_data%d",i), ";hiHF;Events",500,0,5000);
		HFsum_data[i] -> SetMarkerStyle(20);
		HFsum_data[i] -> SetMarkerSize(1.0);
		HFsum_data[i] -> SetMarkerColor(i+1);
		HFsum_data[i] -> SetLabelSize(0.03);
	}

	c_temp -> cd();
	datat_evt -> Draw("hiHF >>+ HFsum_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
	HFsum_data[0] = (TH1D*)gDirectory->Get("HFsum_data0");
	datat_evt -> Draw("hiHF >>+ HFsum_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
	HFsum_data[1] = (TH1D*)gDirectory->Get("HFsum_data1");
	datat_evt -> Draw("hiHF >>+ HFsum_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
	HFsum_data[2] = (TH1D*)gDirectory->Get("HFsum_data2");
	datat_evt -> Draw("hiHF >>+ HFsum_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
	HFsum_data[3] = (TH1D*)gDirectory->Get("HFsum_data3");
	datat_evt -> Draw("hiHF >>+ HFsum_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
	HFsum_data[4] = (TH1D*)gDirectory->Get("HFsum_data4");

	TCanvas *c_HFsum = new TCanvas("c_HFsum", "c_HFsum", 400,400);
	c_HFsum -> SetLogy();
	HFsum_data[0] -> Draw("e1p");
	HFsum_data[1] -> Draw("same&&e1p");
	HFsum_data[2] -> Draw("same&&e1p");
	HFsum_data[3] -> Draw("same&&e1p");
	HFsum_data[4] -> Draw("same&&e1p");

	l1 -> Draw();
	TH1D *HFsum_data_effi[5];
	for(int i=0; i<Ncut; i++){
		//   HFsum_data_effi[i] = new TH1D(Form("HFsum_data_effi%d",i), "", 500,0, 5000);
		HFsum_data_effi[i] = (TH1D*)HFsum_data[i]->Clone(Form("HFsum_data_effi%d",i));
		HFsum_data_effi[i] -> SetTitle(";hiHF;Filter Efficiency");
		if(i!=0)
			HFsum_data_effi[i] -> Divide(HFsum_data[i],HFsum_data[0]);
		HFsum_data_effi[i] -> SetMarkerStyle(20);
		HFsum_data_effi[i] -> SetMarkerSize(1.0);
		HFsum_data_effi[i] -> SetMarkerColor(i+1);
		HFsum_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		HFsum_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
	}

	TCanvas *c_HFsum_effi = new TCanvas("c_HFsum_effi", "c_HFsum_effi", 400, 400);
	HFsum_data_effi[1] -> Draw("elp"); 
	HFsum_data_effi[2] -> Draw("same elp"); 
	HFsum_data_effi[3] -> Draw("same elp"); 
	HFsum_data_effi[4] -> Draw("same elp"); 
	t1 -> Draw();
	l2 -> Draw();

	//======================================
	// hiNpix
	//======================================

	TH1D *HFNpix_data[5];
	for(int i=0; i<Ncut; i++)
	{
		HFNpix_data[i] = new TH1D(Form("HFNpix_data%d",i), ";hiNpix;Events",500,0,10000);
		HFNpix_data[i] -> SetMarkerStyle(20);
		HFNpix_data[i] -> SetMarkerSize(1.0);
		HFNpix_data[i] -> SetMarkerColor(i+1);
	}
	c_temp -> cd();
	datat_evt -> Draw("hiNpix >>+ HFNpix_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
	HFNpix_data[0] = (TH1D*)gDirectory->Get("HFNpix_data0");
	datat_evt -> Draw("hiNpix >>+ HFNpix_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
	HFNpix_data[1] = (TH1D*)gDirectory->Get("HFNpix_data1");
	datat_evt -> Draw("hiNpix >>+ HFNpix_data2", eventCut && "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
	HFNpix_data[2] = (TH1D*)gDirectory->Get("HFNpix_data2");
	datat_evt -> Draw("hiNpix >>+ HFNpix_data3", eventCut && "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
	HFNpix_data[3] = (TH1D*)gDirectory->Get("HFNpix_data3");
	datat_evt -> Draw("hiNpix >>+ HFNpix_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
	HFNpix_data[4] = (TH1D*)gDirectory->Get("HFNpix_data4");

	TCanvas *c_HFNpix = new TCanvas("c_HFNpix", "c_HFNpix", 400,400);
	c_HFNpix -> SetLogy();
	c_HFNpix -> SetLogx();
	HFNpix_data[0] -> Draw("e1p");
	HFNpix_data[1] -> Draw("same&&e1p");
	HFNpix_data[2] -> Draw("same&&e1p");
	HFNpix_data[3] -> Draw("same&&e1p");
	HFNpix_data[4] -> Draw("same&&e1p");
	l1 -> Draw();

	TH1D *HFNpix_data_effi[5];
	for(int i=0; i<Ncut; i++){
		//   HFNpix_data_effi[i] = new TH1D(Form("HFNpix_data_effi%d",i), "", 500,0, 5000);
		HFNpix_data_effi[i] = (TH1D*)HFNpix_data[i]->Clone(Form("HFNpix_data_effi%d",i));
		HFNpix_data_effi[i] -> SetTitle(";hiNpix;Filter Efficiency");
		if(i!=0)
			HFNpix_data_effi[i] -> Divide(HFNpix_data[i],HFNpix_data[0]);
		HFNpix_data_effi[i] -> SetMarkerStyle(20);
		HFNpix_data_effi[i] -> SetMarkerSize(1.0);
		HFNpix_data_effi[i] -> SetMarkerColor(i+1);
		HFNpix_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		HFNpix_data_effi[i] -> SetAxisRange(0.0,1000.0,"X");
	}

	TCanvas *c_HFNpix_effi = new TCanvas("c_HFNpix_effi", "c_HFNpix_effi", 400, 400);
	HFNpix_data_effi[1] -> Draw("elp"); 
	HFNpix_data_effi[2] -> Draw("same elp"); 
	HFNpix_data_effi[3] -> Draw("same elp"); 
	HFNpix_data_effi[4] -> Draw("same elp"); 
	l2 -> Draw();
	t1 -> Draw();



	// ======================================
	// hiZDC
	// ======================================

	TH1D *hiZDC_data[5];
	for(int i=0; i<Ncut; i++){
		hiZDC_data[i] = new TH1D(Form("hiZDC_data%d",i), ";hiZDC;Events",100,0,80000);
		hiZDC_data[i] -> SetMarkerStyle(20);
		hiZDC_data[i] -> SetMarkerSize(0.8);
		hiZDC_data[i] -> SetMarkerColor(i+1);
		hiZDC_data[i] -> SetLabelSize(0.03);
	}

	c_temp -> cd();
	datat_evt -> Draw("hiZDC>>+ hiZDC_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
	hiZDC_data[0] = (TH1D*)gDirectory->Get("hiZDC_data0");
	datat_evt -> Draw("hiZDC>>+ hiZDC_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
	hiZDC_data[1] = (TH1D*)gDirectory->Get("hiZDC_data1");
	datat_evt -> Draw("hiZDC>>+ hiZDC_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
	hiZDC_data[2] = (TH1D*)gDirectory->Get("hiZDC_data2");
	datat_evt -> Draw("hiZDC>>+ hiZDC_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
	hiZDC_data[3] = (TH1D*)gDirectory->Get("hiZDC_data3");
	datat_evt -> Draw("hiZDC>>+ hiZDC_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
	hiZDC_data[4] = (TH1D*)gDirectory->Get("hiZDC_data4");

	TCanvas *c_hiZDC = new TCanvas("c_hiZDC", "c_hiZDC", 400,400);
	c_hiZDC -> SetLogy();
	hiZDC_data[0] -> Draw("e1p");
	hiZDC_data[1] -> Draw("same&&e1p");
	hiZDC_data[2] -> Draw("same&&e1p");
	hiZDC_data[3] -> Draw("same&&e1p");
	hiZDC_data[4] -> Draw("same&&e1p");
	l1 -> Draw();

	TH1D *hiZDC_data_effi[5];
	for(int i=0; i<Ncut; i++){
		hiZDC_data_effi[i] = (TH1D*)hiZDC_data[i]->Clone(Form("hiZDC_data_effi%d",i));
		hiZDC_data_effi[i] -> SetTitle(";hiZDC;Filter Efficiency");
		if(i!=0)
			hiZDC_data_effi[i] -> Divide(hiZDC_data[i],hiZDC_data[0]);
		hiZDC_data_effi[i] -> SetMarkerStyle(20);
		hiZDC_data_effi[i] -> SetMarkerSize(0.8);
		hiZDC_data_effi[i] -> SetMarkerColor(i+1);
		hiZDC_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		hiZDC_data_effi[i] -> SetAxisRange(0.0,80000.0,"X");
	}

	TCanvas *c_hiZDC_effi = new TCanvas("c_hiZDC_effi", "c_hiZDC_effi", 400, 400);
	hiZDC_data_effi[1] -> Draw("elp"); 
	hiZDC_data_effi[2] -> Draw("same elp"); 
	hiZDC_data_effi[3] -> Draw("same elp"); 
	hiZDC_data_effi[4] -> Draw("same elp"); 
	t1 -> Draw();
	l1 -> Draw();



	// ======================================
	// hiEE
	// ======================================

	TH1D *hiEE_data[5];
	for(int i=0; i<Ncut; i++){
		hiEE_data[i] = new TH1D(Form("hiEE_data%d",i), ";hiEE;Events",200,0,4000);
		hiEE_data[i] -> SetMarkerStyle(20);
		hiEE_data[i] -> SetMarkerSize(0.8);
		hiEE_data[i] -> SetMarkerColor(i+1);
		hiEE_data[i] -> SetLabelSize(0.03);
	}

	c_temp -> cd();
	datat_evt -> Draw("hiEE>>+ hiEE_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
	hiEE_data[0] = (TH1D*)gDirectory->Get("hiEE_data0");
	datat_evt -> Draw("hiEE>>+ hiEE_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
	hiEE_data[1] = (TH1D*)gDirectory->Get("hiEE_data1");
	datat_evt -> Draw("hiEE>>+ hiEE_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
	hiEE_data[2] = (TH1D*)gDirectory->Get("hiEE_data2");
	datat_evt -> Draw("hiEE>>+ hiEE_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
	hiEE_data[3] = (TH1D*)gDirectory->Get("hiEE_data3");
	datat_evt -> Draw("hiEE>>+ hiEE_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
	hiEE_data[4] = (TH1D*)gDirectory->Get("hiEE_data4");

	TCanvas *c_hiEE = new TCanvas("c_hiEE", "c_hiEE", 400,400);
	c_hiEE -> SetLogy();
	hiEE_data[0] -> Draw("e1p");
	hiEE_data[1] -> Draw("same&&e1p");
	hiEE_data[2] -> Draw("same&&e1p");
	hiEE_data[3] -> Draw("same&&e1p");
	hiEE_data[4] -> Draw("same&&e1p");
	l1 -> Draw();

	TH1D *hiEE_data_effi[5];
	for(int i=0; i<Ncut; i++){
		hiEE_data_effi[i] = (TH1D*)hiEE_data[i]->Clone(Form("hiEE_data_effi%d",i));
		hiEE_data_effi[i] -> SetTitle(";hiEE;Filter Efficiency");
		if(i!=0)
			hiEE_data_effi[i] -> Divide(hiEE_data[i],hiEE_data[0]);
		hiEE_data_effi[i] -> SetMarkerStyle(20);
		hiEE_data_effi[i] -> SetMarkerSize(0.8);
		hiEE_data_effi[i] -> SetMarkerColor(i+1);
		hiEE_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		hiEE_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
	}

	TCanvas *c_hiEE_effi = new TCanvas("c_hiEE_effi", "c_hiEE_effi", 400, 400);
	hiEE_data_effi[1] -> Draw("elp"); 
	hiEE_data_effi[2] -> Draw("same elp"); 
	hiEE_data_effi[3] -> Draw("same elp"); 
	hiEE_data_effi[4] -> Draw("same elp"); 
	t1 -> Draw();
	l1 -> Draw();


	// ======================================
	// hiET
	// ======================================

	TH1D *hiET_data[5];
	for(int i=0; i<Ncut; i++){
		hiET_data[i] = new TH1D(Form("hiET_data%d",i), ";hiET;Events",100,0,2000);
		hiET_data[i] -> SetMarkerStyle(20);
		hiET_data[i] -> SetMarkerSize(0.8);
		hiET_data[i] -> SetMarkerColor(i+1);
		hiET_data[i] -> SetLabelSize(0.03);
	}

	c_temp -> cd();
	datat_evt -> Draw("hiET>>+ hiET_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
	hiET_data[0] = (TH1D*)gDirectory->Get("hiET_data0");
	datat_evt -> Draw("hiET>>+ hiET_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
	hiET_data[1] = (TH1D*)gDirectory->Get("hiET_data1");
	datat_evt -> Draw("hiET>>+ hiET_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
	hiET_data[2] = (TH1D*)gDirectory->Get("hiET_data2");
	datat_evt -> Draw("hiET>>+ hiET_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
	hiET_data[3] = (TH1D*)gDirectory->Get("hiET_data3");
	datat_evt -> Draw("hiET>>+ hiET_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
	hiET_data[4] = (TH1D*)gDirectory->Get("hiET_data4");

	TCanvas *c_hiET = new TCanvas("c_hiET", "c_hiET", 400,400);
	c_hiET -> SetLogy();
	hiET_data[0] -> Draw("e1p");
	hiET_data[1] -> Draw("same&&e1p");
	hiET_data[2] -> Draw("same&&e1p");
	hiET_data[3] -> Draw("same&&e1p");
	hiET_data[4] -> Draw("same&&e1p");
	l1 -> Draw();

	TH1D *hiET_data_effi[5];
	for(int i=0; i<Ncut; i++){
		hiET_data_effi[i] = (TH1D*)hiET_data[i]->Clone(Form("hiET_data_effi%d",i));
		hiET_data_effi[i] -> SetTitle(";hiET;Filter Efficiency");
		if(i!=0)
			hiET_data_effi[i] -> Divide(hiET_data[i],hiET_data[0]);
		hiET_data_effi[i] -> SetMarkerStyle(20);
		hiET_data_effi[i] -> SetMarkerSize(0.8);
		hiET_data_effi[i] -> SetMarkerColor(i+1);
		hiET_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		hiET_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
	}

	TCanvas *c_hiET_effi = new TCanvas("c_hiET_effi", "c_hiET_effi", 400, 400);
	hiET_data_effi[1] -> Draw("elp"); 
	hiET_data_effi[2] -> Draw("same elp"); 
	hiET_data_effi[3] -> Draw("same elp"); 
	hiET_data_effi[4] -> Draw("same elp"); 
	t1 -> Draw();
	l1 -> Draw();



	// ======================================
	// hiEB
	// ======================================

	TH1D *hiEB_data[5];
	for(int i=0; i<Ncut; i++){
		hiEB_data[i] = new TH1D(Form("hiEB_data%d",i), ";hiEB;Events",200,0,5000);
		hiEB_data[i] -> SetMarkerStyle(20);
		hiEB_data[i] -> SetMarkerSize(0.8);
		hiEB_data[i] -> SetMarkerColor(i+1);
		hiEB_data[i] -> SetLabelSize(0.03);
	}

	c_temp -> cd();
	datat_evt -> Draw("hiEB>>+ hiEB_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
	hiEB_data[0] = (TH1D*)gDirectory->Get("hiEB_data0");
	datat_evt -> Draw("hiEB>>+ hiEB_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
	hiEB_data[1] = (TH1D*)gDirectory->Get("hiEB_data1");
	datat_evt -> Draw("hiEB>>+ hiEB_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
	hiEB_data[2] = (TH1D*)gDirectory->Get("hiEB_data2");
	datat_evt -> Draw("hiEB>>+ hiEB_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
	hiEB_data[3] = (TH1D*)gDirectory->Get("hiEB_data3");
	datat_evt -> Draw("hiEB>>+ hiEB_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
	hiEB_data[4] = (TH1D*)gDirectory->Get("hiEB_data4");

	TCanvas *c_hiEB = new TCanvas("c_hiEB", "c_hiEB", 400,400);
	c_hiEB -> SetLogy();
	hiEB_data[0] -> Draw("e1p");
	hiEB_data[1] -> Draw("same&&e1p");
	hiEB_data[2] -> Draw("same&&e1p");
	hiEB_data[3] -> Draw("same&&e1p");
	hiEB_data[4] -> Draw("same&&e1p");
	l1 -> Draw();

	TH1D *hiEB_data_effi[5];
	for(int i=0; i<Ncut; i++){
		hiEB_data_effi[i] = (TH1D*)hiEB_data[i]->Clone(Form("hiEB_data_effi%d",i));
		hiEB_data_effi[i] -> SetTitle(";hiEB;Filter Efficiency");
		if(i!=0)
			hiEB_data_effi[i] -> Divide(hiEB_data[i],hiEB_data[0]);
		hiEB_data_effi[i] -> SetMarkerStyle(20);
		hiEB_data_effi[i] -> SetMarkerSize(0.8);
		hiEB_data_effi[i] -> SetMarkerColor(i+1);
		hiEB_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		hiEB_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
	}

	TCanvas *c_hiEB_effi = new TCanvas("c_hiEB_effi", "c_hiEB_effi", 400, 400);
	hiEB_data_effi[1] -> Draw("elp"); 
	hiEB_data_effi[2] -> Draw("same elp"); 
	hiEB_data_effi[3] -> Draw("same elp"); 
	hiEB_data_effi[4] -> Draw("same elp"); 
	t1 -> Draw();
	l1 -> Draw();



	// ======================================
	// hiEvtPlanes
	// ======================================

	TH1D *hiEvtPlanes_data[5];
	for(int i=0; i<Ncut; i++){
		hiEvtPlanes_data[i] = new TH1D(Form("hiEvtPlanes_data%d",i), ";hiEvtPlanes;Events",50,-2,2);
		hiEvtPlanes_data[i] -> SetMarkerStyle(20);
		hiEvtPlanes_data[i] -> SetMarkerSize(0.8);
		hiEvtPlanes_data[i] -> SetMarkerColor(i+1);
		hiEvtPlanes_data[i] -> SetLabelSize(0.03);
	}

	c_temp -> cd();
	datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
	hiEvtPlanes_data[0] = (TH1D*)gDirectory->Get("hiEvtPlanes_data0");
	datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
	hiEvtPlanes_data[1] = (TH1D*)gDirectory->Get("hiEvtPlanes_data1");
	datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
	hiEvtPlanes_data[2] = (TH1D*)gDirectory->Get("hiEvtPlanes_data2");
	datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
	hiEvtPlanes_data[3] = (TH1D*)gDirectory->Get("hiEvtPlanes_data3");
	datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
	hiEvtPlanes_data[4] = (TH1D*)gDirectory->Get("hiEvtPlanes_data4");

	TCanvas *c_hiEvtPlanes = new TCanvas("c_hiEvtPlanes", "c_hiEvtPlanes", 400,400);
	c_hiEvtPlanes -> SetLogy();
	hiEvtPlanes_data[0] -> Draw("e1p");
	hiEvtPlanes_data[1] -> Draw("same&&e1p");
	hiEvtPlanes_data[2] -> Draw("same&&e1p");
	hiEvtPlanes_data[3] -> Draw("same&&e1p");
	hiEvtPlanes_data[4] -> Draw("same&&e1p");
	l1 -> Draw();

	TH1D *hiEvtPlanes_data_effi[5];
	for(int i=0; i<Ncut; i++){
		hiEvtPlanes_data_effi[i] = (TH1D*)hiEvtPlanes_data[i]->Clone(Form("hiEvtPlanes_data_effi%d",i));
		hiEvtPlanes_data_effi[i] -> SetTitle(";hiEvtPlanes;Filter Efficiency");
		if(i!=0)
			hiEvtPlanes_data_effi[i] -> Divide(hiEvtPlanes_data[i],hiEvtPlanes_data[0]);
		hiEvtPlanes_data_effi[i] -> SetMarkerStyle(20);
		hiEvtPlanes_data_effi[i] -> SetMarkerSize(0.8);
		hiEvtPlanes_data_effi[i] -> SetMarkerColor(i+1);
		hiEvtPlanes_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		hiEvtPlanes_data_effi[i] -> SetAxisRange(-2.0,2.0,"X");
	}

	TCanvas *c_hiEvtPlanes_effi = new TCanvas("c_hiEvtPlanes_effi", "c_hiEvtPlanes_effi", 400, 400);
	hiEvtPlanes_data_effi[1] -> Draw("elp"); 
	hiEvtPlanes_data_effi[2] -> Draw("same elp"); 
	hiEvtPlanes_data_effi[3] -> Draw("same elp"); 
	hiEvtPlanes_data_effi[4] -> Draw("same elp"); 
	t1 -> Draw();
	l1 -> Draw();

#endif
} 
