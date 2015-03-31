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

void isoVar_mc_data(){
	TH1::SetDefaultSumw2();
	// gStyle->SetOptFit(0);
	// gStyle -> SetTitleYOffset(2.35);
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.05);

//	const double ptbins[] = {20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
	const double ptbins[] = {15,20,30,40,50,60,70,80,90,100,110,120,130,140,160,180,200,240,280};
	const int nptbins = sizeof(ptbins)/sizeof(double) - 1;
	double AvePtBin[nptbins];

	for(int i=0;i<nptbins;i++){
		AvePtBin[i] = (ptbins[i+1]+ptbins[i])/2.0;
	}

	//###################################
	// to merge different pthat samples
	//###################################

	TString treeName = "yJet";
	multiTreeUtil* yJet_mc = new multiTreeUtil();
	yJet_mc->addFile("/u/user/goyeonju/files/yskimfiles/pA/merged_yskim_HiForest_pPb_MIX_AllQCDPhoton_akPu3PF_AfterResCorr_final_l2l3corrTest.root",treeName,"");	
	yJet_mc->AddFriend("tgj");
	multiTreeUtil* yJet_data = new multiTreeUtil();
	yJet_data->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_pPb_DATA_photon30trig_localJEC_v1.root",treeName,"");	
	yJet_data->AddFriend("tgj");

	//###################################
	// flavor 
	//###################################
#if 0
	TCanvas* c_flavor = new TCanvas("c_flavor","c_flavor",400,800);
	c_flavor->Divide(1,2);
	c_flavor->cd(1);
	TH1D* hpt1 = new TH1D("hpt1",";p_{T}^{RECO};Entries",20,0,200);
	TH1D* hpt0 = (TH1D*)hpt1->Clone("hpt0");
	TH1D* hpt2 = (TH1D*)hpt1->Clone("hpt2");
	TH1D* hpt3 = (TH1D*)hpt1->Clone("hpt3");

	yJet -> Draw2(hpt0, "refPt", Form(" photonEt>40 && genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) "),"ptHatWeight*vtxCentWeight");
	yJet -> Draw2(hpt1, "refPt", Form(" photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && refPartonFlv == 21"),"ptHatWeight*vtxCentWeight");
	yJet -> Draw2(hpt2, "refPt", Form(" photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && abs(refPartonFlv)<21 "),"ptHatWeight*vtxCentWeight");
	yJet -> Draw2(hpt3, "refPt", Form(" photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && refPartonFlv < -200"),"ptHatWeight*vtxCentWeight");

	handsomeTH1(hpt0,1);
	handsomeTH1(hpt1,1);
	handsomeTH1(hpt2,2);
	handsomeTH1(hpt3,4);

	hpt0->GetYaxis()->SetTitleOffset(1.8);

	hpt0->DrawCopy("hist");
	hpt1->DrawCopy("same");
	hpt2->DrawCopy("same");
	hpt3->DrawCopy("same");
	jumSun(30,0,30,7400,2);

	c_flavor->cd(2);
	hpt1->Divide(hpt0);
	hpt2->Divide(hpt0);
	hpt3->Divide(hpt0);
	hpt1->SetAxisRange(0,1,"Y");
	hpt1->SetYTitle("Ratio");
	hpt1->DrawCopy();
	hpt2->DrawCopy("same");
	hpt3->DrawCopy("same");
	jumSun(30,0,30,1,2);
#endif

	//###################################
	// to check ptHat spectrum
	//###################################
#if 0
	TCanvas* c_ptHat = new TCanvas("c_ptHat","c_ptHat",400,400);
	TH1D* hptHat = new TH1D("hptHat",";ptHat (GeV);Entries",100,0,500);
	//yJet -> Draw2(hptHat, "ptHat","ptHat>0");
	yJet_mc -> Draw2(hptHat, "pthat","pthat>0","ptHatWeight*vtxCentWeight");
	handsomeTH1(hptHat,1);
	hptHat->DrawCopy();
#endif

	//######################################################
	// jet pt spectrum to compare between MC and DATA 
	//######################################################

	const int nPhoCut = 3;

        TCut etaCut = "";
        if(isforward) etaCut = "-1.66<=eta && eta<-0.47";
        else etaCut = "-0.47<eta && eta<=0.66";
	//TCut kineCut = "dphi>7*3.1415/8. && photonEta<0.50 && photonEta>-1.44";
	TCut kineCut = "dphi>7*3.1415/8. && abs(photonEta)<1.44";
	TCut phoEtCut[nPhoCut] = {"photonEt>40","photonEt>60","photonEt>80"} ;
	TCut decayCut = "sigmaIetaIeta>0.011 && sigmaIetaIeta<0.017";
	TCut sigmaCut = "sigmaIetaIeta<0.01";
	TCut isoCut = "ecalIso<4.2 && hcalIso<2.2 && trackIso<2.0";
	TCut mcCut = "genIso<5 && abs(genMomId)<=22 && genPhotonEt>0";
	const int nTotCut = 3;
	TString str[nPhoCut] = {"photonEt>40", "photonEt>60", "photonEt>80"};
	TCut totCut_data[nTotCut];
	totCut_data[0] = kineCut && etaCut;
	totCut_data[1] = kineCut && sigmaCut && etaCut;
	totCut_data[2] = kineCut && sigmaCut && isoCut && etaCut;
	TCut totCut_mc[nTotCut];
//	totCut_mc[0] = kineCut && mcCut;
//	totCut_mc[1] = kineCut && sigmaCut && mcCut;
	totCut_mc[0] = kineCut && etaCut;
	totCut_mc[1] = kineCut && sigmaCut && etaCut;
	totCut_mc[2] = kineCut && sigmaCut && isoCut && mcCut && etaCut;

	TH1D* hjetpt_data[nTotCut][nPhoCut];
	TH1D* hjetpt_mc[nTotCut][nPhoCut];
	TH1D* hpho_data[nTotCut][nPhoCut];
	TH1D* hpho_mc[nTotCut][nPhoCut];

	TCanvas* c_temp = new TCanvas("c_temp", "jet p_{T} distribution", 1200, 1000); 
	c_temp->cd();
	for(int i=0;i<nTotCut;i++){
		for(int j=0;j<nPhoCut;j++){
			hjetpt_data[i][j] = new TH1D(Form("hjetpt_data_tot%d_pho%d",i,j), Form("hjetpt_data_tot%d_pho%d;p_{T}^{Jet};1/N^{#gamma} dN/dp_{T}^{Jet}",i,j), nptbins, ptbins);
			hjetpt_mc[i][j] = new TH1D(Form("hjetpt_mc_tot%d_pho%d",i,j),  Form("hjetpt_mc_tot%d_pho%d;p_{T}^{Jet};1/N^{#gamma} dN/dp_{T}^{Jet}",i,j), nptbins, ptbins);
			hpho_data[i][j] = new TH1D(Form("hpho_data_tot%d_pho%d",i,j),  Form("hpho_mc_tot%d_pho%d;p_{T}^{#gamma};1/N^{#gamma} dN/dp_{T}^{Jet}",i,j), nptbins, ptbins);
			hpho_mc[i][j] = new TH1D(Form("hpho_mc_tot%d_pho%d",i,j), Form("hpho_mc_tot%d_pho%d;p_{T}^{#gamma};1/N^{#gamma} dN/dp_{T}^{Jet}",i,j), nptbins, ptbins);
		}	
	}

	for(int i=0;i<nTotCut;i++){
		for(int j=0;j<nPhoCut;j++){
			yJet_data -> Draw2(hjetpt_data[i][j], "pt", totCut_data[i] && phoEtCut[j]); 
			yJet_mc -> Draw2(hjetpt_mc[i][j], "pt", totCut_mc[i] && phoEtCut[j],"ptHatWeight*vtxCentWeight"); 
			yJet_data -> Draw2(hpho_data[i][j], "photonEt", totCut_data[i] && phoEtCut[j]); 
			yJet_mc -> Draw2(hpho_mc[i][j], "photonEt", totCut_mc[i] && phoEtCut[j],"ptHatWeight*vtxCentWeight"); 
		}
	}	

	TCanvas* c_cut = new TCanvas("c_cut", "jet p_{T} distribution", 1200, 1000); 
	c_cut->Divide(3,3);

	TLegend *l1 = new TLegend(0.3365615,0.6445304,0.7577623,0.846736,NULL,"brNDC");
	easyLeg(l1);
	for(int i=0;i<nTotCut;i++){
		for(int j=0;j<nPhoCut;j++){
			//if(j==0)cout <<i<<"th, GetTitle : "<<totCut_mc[i].GetTitle()<<endl;
			c_cut->cd(3*i+j+1);
			//gPad->SetLogy();
		
			hjetpt_data[i][j]->Scale(1./hpho_data[i][j]->Integral(),"width");
			hjetpt_mc[i][j]->Scale(1./hpho_mc[i][j]->Integral(),"width");

			//hist cosmetics
			double range = cleverRange(hjetpt_data[i][j],hjetpt_mc[i][j]);
                        hjetpt_data[i][j]->GetYaxis()->SetRangeUser(0.0,range);
                        hjetpt_mc[i][j]->GetYaxis()->SetRangeUser(0.0,range);
			hjetpt_data[i][j]->SetMarkerStyle(20);
			hjetpt_mc[i][j]->SetLineColor(4);
			hjetpt_mc[i][j]->Draw("hist e");
			hjetpt_data[i][j]->Draw("same p");
			if(i==0 && j==0) {
				l1->AddEntry(hjetpt_data[0][0],"Data","p");
				l1->AddEntry(hjetpt_mc[0][0],"MC","l");
				l1->Draw("same");
			}
			drawText(str[j],0.6,0.75);			
			/*	
				if(j==0){
				TLegend* leg = new TLegend(0.3,0.3,0.9,0.9);
				leg->AddEntry((TObject*)0,Form("%s",totCut_data[i].GetTitle()),"");
				leg->SetTextColor(kBlack);
				leg->Draw();
				}	
				*/
			//if(j==0)drawLongText(0.5,0.5,0.7,0.7,Form("%s",totCut_data[i].GetTitle()),kBlack);
			//if(j==1)drawLongText(0.5,0.5,0.7,0.7,Form("%s",totCut_mc[i].GetTitle()),kBlue);
		}
	}
//	c_cut -> Update();
	c_cut ->SaveAs(Form("pdf/yJetptdist_isforward%d.pdf",isforward));
}

