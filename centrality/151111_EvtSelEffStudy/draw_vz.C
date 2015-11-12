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

void draw_vz(const char* fname="/u/user/goyeonju/files/centrality/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611_CMSSW5320_byYJ.root", TString type="2011PbPb_5320", bool isAOD=0)
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

    TCanvas *c[3];
    TCanvas *c2D[3];
    TH1D* h[3];
    TH2D* h2D[3];
    for(int i=0;i<3;i++){
        c[i] = new TCanvas(Form("c%d",i),"",300,300);
        c2D[i] = new TCanvas(Form("c2D%d",i),"",300,300);
        h[i] = new TH1D(Form("h%d",i),";z_{vtx} (cm);",50,-50,50);
        h2D[i] = new TH2D(Form("h2D%d",i),";x_{vtx} (cm);y_{vtx} (cm)",100,0.0,0.4,100,0.0,0.4);
    }
    t_evt->Draw("vz>>h0","HLT_HIMinBiasHfOrBSC_v1==1");
    t_evt->Draw("vz>>h1","HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    t_evt->Draw("vz>>h2","HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==0"); 
    t_evt->Draw("vy:vx>>h2D0","HLT_HIMinBiasHfOrBSC_v1==1");
    t_evt->Draw("vy:vx>>h2D1","HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    t_evt->Draw("vy:vx>>h2D2","HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==0"); 

    for(int i=0;i<3;i++){
        c[i]->cd();
        h[i]->Draw("p");
        c[i]->SaveAs(Form("pdf/vz_pprimaryVertexFilter%d.pdf",i));

        c2D[i]->cd();
        h2D[i]->Draw("colz");
        gPad->SetLogz();
        c2D[i]->SaveAs(Form("pdf/vyvx_pprimaryVertexFilter%d.pdf",i));
    }
}
