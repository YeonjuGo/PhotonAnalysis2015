// comparison centrality objects between runs 
// Prompt Reco vs. Private Reco vs. Express
// Author : Yeonju Go
// 08 Dec 2015

//private setup
#include "../../yjUtility.h"

//don't use kOrange, kGreen, kYellow!!!!!
const Int_t colhere[10] = {1, kBlue, kRed, kCyan+2, kGray+1, kMagenta+2, kYellow+1,kGreen+3,   kMagenta,  kCyan-9};
//const int colhere[] = {1,2,kViolet,4,kOrange,kGreen,28};
//const int colhere[] = {1,4,9,6,46,kOrange,kViolet,kOrange+10};
//const int colhere[] = {1,2,4,9,6,46,kOrange,kViolet,kOrange+10};
void draw_compare_runs_1D(string var="hiHF");
void compare_runs_1D_3types(){
    TStopwatch timer;
    timer.Start();
    draw_compare_runs_1D("hiHF");
#if 0
    draw_compare_runs_1D("zVtx");//change many things !! 
    draw_compare_runs_1D("xVtx");
    draw_compare_runs_1D("yVtx");
    draw_compare_runs_1D("hiNpix");
    draw_compare_runs_1D("hiHF");
    draw_compare_runs_1D("hiNpix");
    draw_compare_runs_1D("hiNtracks");
    draw_compare_runs_1D("sumPtVtx");
    draw_compare_runs_1D("hiBin");
    draw_compare_runs_1D("hiEE");
    draw_compare_runs_1D("hiEB");
#endif
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
    const char* dir = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/";
    const int Nrun = 7;
    string run[] = {"2811","2694","2695","2816","2620","2656","3233"};
    //const char* trigCut = "HLT_HIL1MinimumBiasHF1AND_v1";
    //TCut trigCut = "HLT_HIL1MinimumBiasHF1AND_v1 && pcollisionEventSelection";
    TCut trigCut = "hiHF>220 && HLT_HIL1MinimumBiasHF1AND_v1 && pprimaryVertexFilter && phfCoincFilter3";
    TCut trigCutFor3233 = "hiHF>220 && HLT_HIL1MinimumBiasHF2AND_part2_v1 && pprimaryVertexFilter && phfCoincFilter3";
    const char* cap= "";

    int nBin = 100;
    if(var=="hiBin") nBin=200;
    double xmax = 6000;
    if(var == "hiNpix") xmax = 50000;
    if(var == "hiNtracks" || var == "sumPtVtx") xmax = 3500;
    if(var == "hiBin") xmax = 200;
    if(var == "hiEE") xmax = 3000;
    if(var == "hiEB") xmax = 4000;
    if(var == "xVtx" || var == "yVtx") xmax = 0.20;

    int nEvents[Nrun];
    TH1D* h[Nrun];
    TLegend* l1 = new TLegend(0.7,0.65,0.95,0.95);
    legStyle(l1);
    for(int i=0; i<Nrun; ++i){
        TFile* f1;
        if(run[i]=="2811" ||run[i]=="2816") f1 = TFile::Open(Form("%sHIMinimumBias2_run26%s.root",dir,run[i].data()));
        else if(run[i]=="2548") f1 = TFile::Open(Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HIForestMinbiasUPC_run262548.root"));
        //else if(run[i]=="2620" || run[i]=="2656") f1 = TFile::Open(Form("/data/velicanu/store/group/phys_heavyions/velicanu/forest/HIRun2015/ReReco/Merged/HIMinimumBias2.26%s.root",run[i].data()));
        else if(run[i]=="3233") f1 = TFile::Open(Form("root://eoscms//eos/cms/store/group/cmst3/group/hintt/mverweij/PbPbReco/Forest/000/263/233/merge/HiForestRun263233PrivateRecoMB.root"));
        else if(run[i]=="2620"||run[i]=="2656") f1 = TFile::Open(Form("%sHIForestExpress_run26%s.root",dir,run[i].data()));
        else f1 = TFile::Open(Form("%sHiForestPromptReco_26%s.root",dir,run[i].data()));
        if(f1->IsZombie()) { cout << run[i].data() << " doesn't exist!! " << endl; continue;}
        cout << "Open file : " << f1->GetName() << endl;
        TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
        TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
        TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
        TTree* t1_track= (TTree*) f1 -> Get("anaTrack/trackTree");
        t1->AddFriend(t1_skim);
        t1->AddFriend(t1_hlt);
        t1->AddFriend(t1_track);
        if(var == "zVtx") h[i] = new TH1D(Form("h%d",i), Form(";%s (cm);Event Fraction",var.data()), nBin, -25, 25); // for zVtx 
        else if(var == "yVtx"|| var=="xVtx") h[i] = new TH1D(Form("h%d",i), Form(";%s (cm);Event Fraction",var.data()), nBin, 0.06, 0.12); // for x,yVtx 
        else h[i] = new TH1D(Form("h%d",i), Form(";%s;Event Fraction",var.data()), nBin, 0, xmax);
        if(run[i]=="3233") {
            t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCutFor3233);
            nEvents[i] = t1->GetEntries(trigCutFor3233);
        }
        else {
            t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCut);
            nEvents[i] = t1->GetEntries(trigCut);
        } 
            cout << run[i].data() << " # of events : " << nEvents[i] << endl;
    }

    ///////////////////////////////////////////////////////
    // vertex study

    if(var == "xVtx" || var == "yVtx" || var == "zVtx"){ 
        TH1D* hMean = new TH1D("hMean", Form(";Run;<%s>",var.data()),Nrun,0,Nrun+1);
        TH1D* hSig = new TH1D("hMean", Form(";Run;#sigma(%s)",var.data()),Nrun,0,Nrun+1);
        hMean->SetMarkerStyle(20);
        hSig->SetMarkerStyle(20);
        for(int i=0; i<Nrun; ++i){
            hMean->SetBinContent(i+1,  h[i]->GetMean());
            hMean->SetBinError(i+1,  h[i]->GetMeanError());
            hSig->SetBinContent(i+1, h[i]->GetRMS());
            hSig->SetBinError(i+1, h[i]->GetRMSError());
        }
        TCanvas* cM = new TCanvas("cM","",500,500);
        hMean->Draw();
        cM->SaveAs(Form("pdf/compareBtwRuns_run26%s_%s_mean_%s.pdf",run[0].data(),var.data(),cap));
        TCanvas* cS = new TCanvas("cS","",500,500);
        hSig->Draw();
        cS->SaveAs(Form("pdf/compareBtwRuns_run26%s_%s_sigma_%s.pdf",run[0].data(),var.data(),cap));
    }


    //////////////////////////////////////////////////////
    // DRAW!!
    TCanvas* c1 = new TCanvas("c2","",500,500);
    if(var != "hiBin") c1->SetLogy();
    for(int i=0; i<Nrun; ++i){
        h[i]->Scale(1./nEvents[i]); // this is for all!! 
        //h[i]->Scale(1./h[i]->Integral("width"));
        SetHistColor(h[i],colhere[i]);
        if(i==0) h[i]->DrawCopy("hist");
        //else if(i>=3) h[i]->DrawCopy("p same");
        else h[i]->DrawCopy("hist same");
        l1->AddEntry(h[i], Form("Run 26%s",run[i].data()),"pl"); 
    }
    l1->Draw();
    c1->SaveAs(Form("pdf/compareBtwRuns_run26%s_%s%s.pdf",run[0].data(),var.data(),cap));

    if(var=="hiBin") h[0]->Rebin(10);
    else h[0]->Rebin(5);
    TCanvas* c2 = new TCanvas("c3","", 500, 500); 
    for(int i=1; i<Nrun; ++i){
        if(var=="hiBin") h[i]->Rebin(10);
        else h[i]->Rebin(5);
        h[i]->Divide(h[0]);
        h[i]->SetMarkerStyle(20);
        h[i]->SetMarkerSize(0.85);
        h[i]->GetYaxis()->SetRangeUser(0.6,1.3);
        h[i]->SetTitle(Form(";%s;Run XXX / Run 26%s",var.data(),run[0].data()));
        if(i==1) h[i]->DrawCopy("pe");
        else h[i]->DrawCopy("pe same");
    }
    if(var == "zVtx") jumSun(-25,1,25,1);//for zVtx
    else jumSun(0,1,xmax,1);
    c2->SaveAs(Form("pdf/compareBtwRuns_ratio_run26%s_%s%s.pdf",run[0].data(),var.data(),cap));

   

} 

