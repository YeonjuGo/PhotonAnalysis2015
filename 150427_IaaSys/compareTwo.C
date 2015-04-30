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

void compareTwo(int icoll=1){
	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.05);

    TString coll_str="";
    if(icoll==0) coll_str="PbPb";
    else if(icoll==1) coll_str="pp13099";
    else if(icoll==2) coll_str="pp10030";
    //icoll 0=PbPb, 1=pp 30-100%, 2=pp 0-30%

    TLegend *l1 = new TLegend(0.52, 0.55, 0.86, 0.92);
    if(icoll==0) l1->AddEntry((TObject*)0, "PbPb", "");
    else if(icoll==1) l1->AddEntry((TObject*)0, "pp 30-100 %", "");
    else if(icoll==2) l1->AddEntry((TObject*)0, "pp 0-30 %", "");

    TFile * f1;//old
    TFile * f2;//new
    if(icoll==0){//PbPb
        f1 = new TFile("/home/goyeonju/CMS/2015/gammaJetSystematics/resultHistograms/resultHistograms_nominal_vtxCentWeighted.root");
        f2 = new TFile("/home/goyeonju/CMS/2015/gammaJetSystematics/resultHistograms_yj/resultHistograms_nominal_vtxCentWeighted.root");
    } else if(icoll==1){//pp 30-100 % 
        f1 = new TFile("/home/goyeonju/CMS/2015/gammaJetSystematics/resultHistograms/resultHistograms_ppSmeared13099.root");
        f2 = new TFile("/home/goyeonju/CMS/2015/gammaJetSystematics/resultHistograms_yj/resultHistograms_ppSmeared13099.root"); 
    } else if(icoll==2){//pp 0-30 % 
        f1 = new TFile("/home/goyeonju/CMS/2015/gammaJetSystematics/resultHistograms/resultHistograms_ppSmeared10030.root");
        f2 = new TFile("/home/goyeonju/CMS/2015/gammaJetSystematics/resultHistograms_yj/resultHistograms_ppSmeared10030.root");
    }

//////////////////////////////////////////////////// 
//////////////////////////////////////////////////// 
// dNdJetPt_pp_ptBin1

    TH1D* h_IaaBin_pp[4][2];
    TH1D* h_IaaBin_pp_ratio[4];
    for(int i=0; i<4; i++){
        h_IaaBin_pp[i][0] = (TH1D*) f1 -> Get(Form("dNdJetPt_IaaBin_pp_ptBin%d",i+1));
        h_IaaBin_pp[i][1] = (TH1D*) f2 -> Get(Form("dNdJetPt_IaaBin_pp_ptBin%d",i+1));
        h_IaaBin_pp_ratio[i] = (TH1D*) h_IaaBin_pp[i][0]->Clone(Form("IaaBin_pp%d_ratio",i+1));
        h_IaaBin_pp_ratio[i]->Clear();
        h_IaaBin_pp_ratio[i]->Divide(h_IaaBin_pp[i][1],h_IaaBin_pp[i][0]);
        h_IaaBin_pp_ratio[i]->SetTitle(";Jet p_{T} (GeV);Ratio (New/Old)");
    }

    //draw
	TCanvas* c1=new TCanvas("c1","",1300,500);
    c1->Divide(4,2);
    for(int i=0; i<4; i++){
        c1->cd(i+1);
        h_IaaBin_pp[i][1]->SetMarkerStyle(20);
        h_IaaBin_pp[i][1]->Draw("p");
        h_IaaBin_pp[i][0]->Draw("same hist");
        c1->cd(i+5);
        h_IaaBin_pp_ratio[i]->SetMarkerStyle(20);
        h_IaaBin_pp_ratio[i]->Draw();
        jumSun(30,1,200,1,4);
    }
    c1->cd(1);
    l1->AddEntry(h_IaaBin_pp[0][0],"Old","l");
    l1->AddEntry(h_IaaBin_pp[0][1],"New","p");
    l1->Draw();
    c1 ->SaveAs(Form("pdf/compare_dNdJetPtIaaBin_pp_%s.pdf",coll_str.Data()));


//////////////////////////////////////////////////// 
//////////////////////////////////////////////////// 
// dNdJetPt_IaaBin_pbpb_centralityBin1_ptBin1 

    TH1D* h_IaaBin_pbpb_centralityBin1[4][2];
    TH1D* h_IaaBin_pbpb_centralityBin1_ratio[4];
    for(int i=0; i<4; i++){
        h_IaaBin_pbpb_centralityBin1[i][0] = (TH1D*) f1 -> Get(Form("dNdJetPt_IaaBin_pbpb_centralityBin1_ptBin%d",i+1));
        h_IaaBin_pbpb_centralityBin1[i][1] = (TH1D*) f2 -> Get(Form("dNdJetPt_IaaBin_pbpb_centralityBin1_ptBin%d",i+1));
        h_IaaBin_pbpb_centralityBin1_ratio[i] = (TH1D*) h_IaaBin_pbpb_centralityBin1[i][0]->Clone(Form("IaaBin_pbpb_centralityBin1%d_ratio",i+1));
        h_IaaBin_pbpb_centralityBin1_ratio[i]->Clear();
        h_IaaBin_pbpb_centralityBin1_ratio[i]->Divide(h_IaaBin_pbpb_centralityBin1[i][1],h_IaaBin_pbpb_centralityBin1[i][0]);
        h_IaaBin_pbpb_centralityBin1_ratio[i]->SetTitle(";Jet p_{T} (GeV);Ratio (New/Old)");
    }

    //draw
	TCanvas* c2=new TCanvas("c2","",1300,500);
    c2->Divide(4,2);
    for(int i=0; i<4; i++){
        c2->cd(i+1);
        h_IaaBin_pbpb_centralityBin1[i][1]->SetMarkerStyle(20);
        h_IaaBin_pbpb_centralityBin1[i][1]->Draw("p");
        h_IaaBin_pbpb_centralityBin1[i][0]->Draw("same hist");
        c2->cd(i+5);
        h_IaaBin_pbpb_centralityBin1_ratio[i]->SetMarkerStyle(20);
        h_IaaBin_pbpb_centralityBin1_ratio[i]->Draw();
        jumSun(30,1,200,1,4);
    }
    c2->cd(1);
    l1->Draw();
    c2 ->SaveAs(Form("pdf/compare_dNdJetPtIaaBin_pbpb_centralityBin1_%s.pdf",coll_str.Data()));


//////////////////////////////////////////////////// 
//////////////////////////////////////////////////// 
// dNdJetPt_IaaBin_pbpb_centralityBin2_ptBin1 

    TH1D* h_IaaBin_pbpb_centralityBin2[4][2];
    TH1D* h_IaaBin_pbpb_centralityBin2_ratio[4];
    for(int i=0; i<4; i++){
        h_IaaBin_pbpb_centralityBin2[i][0] = (TH1D*) f1 -> Get(Form("dNdJetPt_IaaBin_pbpb_centralityBin2_ptBin%d",i+1));
        h_IaaBin_pbpb_centralityBin2[i][1] = (TH1D*) f2 -> Get(Form("dNdJetPt_IaaBin_pbpb_centralityBin2_ptBin%d",i+1));
        h_IaaBin_pbpb_centralityBin2_ratio[i] = (TH1D*) h_IaaBin_pbpb_centralityBin2[i][0]->Clone(Form("IaaBin_pbpb_centralityBin2%d_ratio",i+1));
        h_IaaBin_pbpb_centralityBin2_ratio[i]->Clear();
        h_IaaBin_pbpb_centralityBin2_ratio[i]->Divide(h_IaaBin_pbpb_centralityBin2[i][1],h_IaaBin_pbpb_centralityBin2[i][0]);
        h_IaaBin_pbpb_centralityBin2_ratio[i]->SetTitle(";Jet p_{T} (GeV);Ratio (New/Old)");
    }

    //draw
	TCanvas* c3=new TCanvas("c3","",1300,500);
    c3->Divide(4,2);
    for(int i=0; i<4; i++){
        c3->cd(i+1);
        h_IaaBin_pbpb_centralityBin2[i][1]->SetMarkerStyle(20);
        h_IaaBin_pbpb_centralityBin2[i][1]->Draw("p");
        h_IaaBin_pbpb_centralityBin2[i][0]->Draw("same hist");
        c3->cd(i+5);
        h_IaaBin_pbpb_centralityBin2_ratio[i]->SetMarkerStyle(20);
        h_IaaBin_pbpb_centralityBin2_ratio[i]->Draw();
        jumSun(30,1,200,1,4);
    }
    c3->cd(1);
    l1->Draw();
    c3 ->SaveAs(Form("pdf/compare_dNdJetPtIaaBin_pbpb_centralityBin2_%s.pdf",coll_str.Data()));
}

