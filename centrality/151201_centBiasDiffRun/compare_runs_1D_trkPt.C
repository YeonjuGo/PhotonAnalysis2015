// comparison centrality objects between runs 
// Author : Yeonju Go
// 30 Nov 2015

//basic c++ header, string ...
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
#include <TCut.h>
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
#include <string>
//private setup
#include "yjUtility.h"

const int colhere[] = {1,2,4,9,6,46,kOrange};
void draw_compare_runs_1D(string var="hiHF");
void compare_runs_1D_trkPt(){
    TStopwatch timer;
    timer.Start();
    draw_compare_runs_1D("trkPt");
    //draw_compare_runs_1D("zVtx");//change many things !! 

    timer.Stop();
    cout<<"Macro finished: "<<endl;
    cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
    cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
}

void draw_compare_runs_1D(string var){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2(); 
 
    ///////////////////// for larger sample
    const char* dir = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIMinimumBias2/Merged/";
    const int Nrun = 6;
    string run[] = {"695","548","620","656","694","816"};
    //string run[] = {"695","548","620","656","694","811","816"};
    //const char* trigCut = "HLT_HIL1MinimumBiasHF1AND_v1";
    TCut trigCut = "HLT_HIL1MinimumBiasHF1AND_v1 && pprimaryVertexFilter && phfCoincFilter3";
    //const char* cap= "_yLinear";

    double trkPtBin[] = {0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0,12.0,16.0,20.0};
    int NtrkPtBin = sizeof(trkPtBin)/sizeof(double)-1;

    int nBin = 100;
    double xmax = 40;

    int nEvents[Nrun];

    TH1D* h[Nrun];
    TLegend* l1 = new TLegend(0.5,0.6,0.95,0.9);
    legStyle(l1);
    for(int i=0; i<Nrun; ++i){
        TFile* f1;
        if(run[i]=="811") f1 = TFile::Open(Form("/afs/cern.ch/work/y/ygo/public/PFphoton/HIForestMinimumBias_run262811_file1.root"));
        else if(run[i]=="816") f1 = TFile::Open(Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIMinimumBias_run262816.root"));
        else if(run[i]=="548") f1 = TFile::Open(Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HIForestMinbiasUPC_run262548.root"));
        else f1 = TFile::Open(Form("%sHIForestExpress_run262%s.root",dir,run[i].data()));
        if(f1->IsZombie()) { cout << run[i].data() << " doesn't exist!! " << endl; continue;}
        cout << "Open file : " << f1->GetName() << endl;
        TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
        TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
        TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
        TTree* t1_track= (TTree*) f1 -> Get("anaTrack/trackTree");
        t1->AddFriend(t1_skim);
        t1->AddFriend(t1_hlt);
        t1->AddFriend(t1_track);
        nEvents[i] = t1->GetEntries(trigCut);
        cout << run[i].data() << " # of events : " << nEvents[i] << endl;
        h[i] = new TH1D(Form("h%d",i), Form(";%s (GeV);1/N_{evt}dN/dp_{T}",var.data()), NtrkPtBin, trkPtBin);
        
        t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCut && "highPurity==1");//for tracks
    }

    TCanvas* c1 = new TCanvas("c2","",500,500);
    c1->SetLogy();
    for(int i=0; i<Nrun; ++i){
        h[i]->Scale(1./nEvents[i],"width"); // this is for all!! 
        SetHistColor(h[i],colhere[i]);
        if(i==0) h[i]->DrawCopy("hist");
        else h[i]->DrawCopy("hist same");
        l1->AddEntry(h[i], Form("Run 262%s",run[i].data()),"pl"); 
    }
    l1->Draw();
    c1->SaveAs(Form("pdf/compareBtwRuns_run262%s_%s%s.pdf",run[0].data(),var.data(),cap));

    //h[0]->Rebin(4);
    TCanvas* c2 = new TCanvas("c3","", 500, 500); 
    for(int i=1; i<Nrun; ++i){
        //h[i]->Rebin(4);
        h[i]->Divide(h[0]);
        h[i]->SetMarkerStyle(20);
        h[i]->SetMarkerSize(0.85);
        h[i]->GetYaxis()->SetRangeUser(0.6,1.3);
        h[i]->SetTitle(Form(";%s;Run XXX / Run 262695",var.data()));
        if(i==1) h[i]->DrawCopy("pe");
        else h[i]->DrawCopy("pe same");
    }
    if(var == "zVtx") jumSun(-25,1,25,1);//for zVtx
    else jumSun(0,1,xmax,1);
    c2->SaveAs(Form("pdf/compareBtwRuns_ratio_run262%s_%s%s.pdf",run[0].data(),var.data(),cap));
} 

