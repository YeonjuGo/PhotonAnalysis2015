// Author Yeonju Go
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
#include "../../yjUtility.h"

const double dy= 0.5;
const int Ncut = 5;
void Get1DEffPlots(TTree* t_evt=0, TString v1="hiHF",int xbin=200, double xmin=0, double xmax=4500, TCut cut="", TCanvas* c_tot=0, TString cap="", bool isPassed=0, double eff_ymin=0.50, bool isAOD=0);

void draw_vz(const char* fname="/u/user/goyeonju/files/centrality/Centrality_officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root", TString type="HYDJET_5320", bool isAOD=0)
{
    const TCut runCut = "run==181611";
    const TCut lumiCut = "lumi>=1 && lumi<=895";
    const TCut eventCut = runCut && lumiCut;
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TFile *fin = new TFile(fname);
    TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
    t_evt -> AddFriend(t_hlt);
    t_evt -> AddFriend(t_skim);
    double Nevt_t = t_evt -> GetEntries();
    cout << "# of DATA events = " << Nevt_t << endl;


    TCanvas *c_tot = new TCanvas("c_tot", "c_tot", 900,300);
    c_tot->Divide(3,1);
    TH1D* h1 = new TH1D("h1",";vz;",50,-50,50);
    TH1D* h2 = new TH1D("h2",";vz;",50,-50,50);
    TH1D* h3 = new TH1D("h3",";vz;",50,-50,50);
    t_evt->Draw("vz>>h1","HLT_HIMinBiasHfOrBSC_v1");
    t_evt->Draw("vz>>h2","HLT_HIMinBiasHfOrBSC_v1 && pprimaryVertexFilter==1");
    t_evt->Draw("vz>>h3","HLT_HIMinBiasHfOrBSC_v1 && pprimaryVertexFilter==0");
    c_tot->cd(1);
    h1->Draw();
    c_tot->cd(2);
    h2->Draw();
    c_tot->cd(3);
    h3->Draw();
}

