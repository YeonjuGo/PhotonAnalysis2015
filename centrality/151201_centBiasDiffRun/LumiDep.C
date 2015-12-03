// to draw the distribution of centrality variables
// with the different event selection filters in one canvas.
// should modify evtfilter definition part
//
// Author Yeonju Go
// 23 Nov 2015
//
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "stdio.h"
#include <iostream>
#include "yjUtility.h"

void GetLumiDepPlots(string run="2620", string var="hiBin", int nBin_var=10, double var_i=0, double var_f=10, int nBin_lumi=120, double lumi_i=100, double lumi_f=220, TCut cut="", const char* cap="");
void LumiDep()
{
   // GetLumiDepPlots("2620", "hiBin",5,0,5, 120,100,220, "", "");
   // GetLumiDepPlots("2656", "hiBin",5,0,5, 120,1,178, "", "");
    GetLumiDepPlots("2694", "hiBin",5,0,5, 120,71,229, "", "");
    GetLumiDepPlots("2695", "hiBin",5,0,5, 120,1,227, "", "");
   // GetLumiDepPlots("2816", "hiBin",5,0,5, 120,1,449, "", "");
}

void GetLumiDepPlots(string run, string var, int nBin_var, double var_i, double var_f, int nBin_lumi, double lumi_i, double lumi_f, TCut cut, const char* cap)
{
    TCut trig= "HLT_HIL1MinimumBiasHF1AND_v1 && pprimaryVertexFilter && phfCoincFilter3";
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle -> SetOptStat(1);
    SetHistTitleStyle(0.06,0.04);
    SetyjPadStyle();

    const char* fname="";
    if(run=="2810" || run=="2816"){
        fname = Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIMinimumBias2_run26%s.root", run.data());
    } else if(run=="2620" || run=="2656" || run=="2694" || run=="2695"){
        fname = Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIMinimumBias2/Merged/HIForestExpress_run26%s.root", run.data());
    } else if(run=="2548"){
        fname = Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HIForestMinbiasUPC_run262548.root");
    }
    TFile* fin = TFile::Open(fname);
    TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
    t_evt -> AddFriend(t_hlt);
    t_evt -> AddFriend(t_skim);

    double Nevt_t = t_evt -> GetEntries(trig);
    cout << "Trigger is '" << trig.GetTitle() << "' and # of events = " << Nevt_t << endl;

    TCanvas *c1 = new TCanvas(Form("c1_%s_%s_run26%s",var.data(),cap,run.data()), "c_tot", 600,500);
    TH2D* h2D = new TH2D(Form("h2D_%s_lumi",var.data()),Form(";Lumi;%s", var.data()), nBin_lumi, lumi_i, lumi_f, nBin_var, var_i, var_f);
    TProfile* prof;
    t_evt->Draw(Form("%s:lumi>>+%s",var.data(),h2D->GetName()), trig && cut && Form("%s>=%f && %s<=%f",var.data(),(float)var_i,var.data(),(float)var_f), "colz");
    h2D=(TH2D*)gDirectory->Get(h2D->GetName());
    //gPad->SetLogz();
    prof = h2D->ProfileX();
    prof->SetMarkerStyle(20);
    c1->cd();
    h2D->Draw("colz");
    prof->Draw("same");
    drawText(Form("run26%s",run.data()),0.3,0.95);
    c1->SaveAs(Form("pdf/lumiDep_2D_%s_%s_run26%s.pdf",var.data(),cap,run.data()));
}
