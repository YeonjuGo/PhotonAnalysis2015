// Author : Yeonju Go
#include "../../yjUtility.h"

const Int_t colhere[10] = {1, kBlue, kRed, kCyan+2, kGray+1, kMagenta+2, kYellow+1,kGreen+3, kMagenta,  kCyan-9};
//const int colhere[] = {1,4,9,2,6,46,28,kViolet,kOrange+10};
void draw_compare_runs_1D(string var="hiHF");
void compare_runs_1D_express(){
    TStopwatch timer;
    timer.Start();
    draw_compare_runs_1D("hiHF");
#if 0
    draw_compare_runs_1D("hiHF");
    draw_compare_runs_1D("hiBin");
    draw_compare_runs_1D("hiNpix");
    draw_compare_runs_1D("hiNtracks");
    draw_compare_runs_1D("zVtx");//change many things !! 
    draw_compare_runs_1D("hiEE");
    draw_compare_runs_1D("hiEB");
    //draw_compare_runs_1D("sumPtVtx");
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
#if 1 
    ////////////////////// for express
    const char* dir = "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/";
    const int Nrun = 4;
    string run[] = {"2811","3035","3233","3234"};
    //string run[] = {"2811","3035","3233","3234","3284","3293","3349","3362"};
    //string run[] = {"2811","3035","3233","3234","3284","3293","3349","3362"};
    //string run[] = {"2811","3284","3293","3333","3349","3355","3358","3362"};
    //string run[] = {"3284","3286","3293", "3322","3333","3349","3355","3358","3362"};
    //string run[] = {"3284","3286","3293", "3322","3333","3349","3355","3356","3357","3358","3359","3362"};
    //string run[] = {"2811","2816", "2817", "2818", "2819"};
    //string run[] = {"2811", "2815", "2816", "2817", "2818", "2819"};
    //string run[] = {"2768","2620", "2656", "2695", "2726", "2811", "2816", "2837","3035"};
    const char* trigCutHF1 = "(HLT_HIL1MinimumBiasHF1ANDExpress_v1) && pcollisionEventSelection"; 
    const char* trigCutHF2 = "(HLT_HIL1MinimumBiasHF1ANDExpress_v1 || HLT_HIL1MinimumBiasHF2ANDExpress_v1) && pcollisionEventSelection"; 
    //const char* trigCut = "HLT_HIL1MinimumBiasHF1ANDExpress_v1 && pprimaryVertexFilter && phfCoincFilter3";
    const char* cap= "";

#endif
    int nBin = 50;
    if(var=="hiBin") nBin=50;
    double xmax = 6000;
    if(var == "hiNpix") xmax = 50000;
    if(var == "hiNtracks" || var == "sumPtVtx") xmax = 3500;
    if(var == "hiBin") xmax = 200;
    if(var == "hiEE") xmax = 3000;
    if(var == "hiEB") xmax = 4000;

    int nEvents[Nrun];
    TH1D* h[Nrun];
    TLegend* l1 = new TLegend(0.7,0.65,0.95,0.95);
    legStyle(l1);
    for(int i=0; i<Nrun; ++i){
        TFile* f1;
        if(run[i]=="2837") f1 = TFile::Open(Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262837_preHBHE.root"));
        else if(run[i]=="2620") f1 = TFile::Open(Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262620-v6.root"));
        else f1 = TFile::Open(Form("%sHIForestExpress_run26%s.root",dir,run[i].data()));
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
        else h[i] = new TH1D(Form("h%d",i), Form(";%s;Event Fraction",var.data()), nBin, 0, xmax);
        if(run[i]=="2811") {
            nEvents[i] = t1->GetEntries(trigCutHF1);
            t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCutHF1);
        } else { 
            nEvents[i] = t1->GetEntries(trigCutHF2);
            t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCutHF2);
        }
        cout << run[i].data() << " # of events : " << nEvents[i] << endl;
    }
    TCanvas* c1 = new TCanvas("c2","",500,500);
    if(var != "hiBin") c1->SetLogy();
    for(int i=0; i<Nrun; ++i){
        h[i]->Scale(1./nEvents[i]); // this is for all!! 
        //h[i]->Scale(1./h[i]->Integral("width"));
        SetHistColor(h[i],colhere[i]);
        if(i==0) h[i]->DrawCopy("hist");
        else h[i]->DrawCopy("hist same");
        l1->AddEntry(h[i], Form("Run 26%s",run[i].data()),"pl"); 
    }
    cout << "ss" << endl;
    l1->Draw();
    c1->SaveAs(Form("pdf/compareBtwRuns_refrun26%s_%s%s.pdf",run[0].data(),var.data(),cap));

    h[0]->Rebin(5);
    TCanvas* c2 = new TCanvas("c3","", 500, 500); 
    cout << "ss" << endl;
    for(int i=1; i<Nrun; ++i){
        h[i]->Rebin(5);
        h[i]->Divide(h[0]);
        h[i]->SetMarkerStyle(20);
        h[i]->SetMarkerSize(0.85);
        h[i]->GetYaxis()->SetRangeUser(0.5,1.3);
        h[i]->SetTitle(Form(";%s;Run XXX / Run 26%s",var.data(),run[0].data()));
        if(i==1) h[i]->DrawCopy("pe");
        else h[i]->DrawCopy("pe same");
    }
    cout << "ss" << endl;
    if(var == "zVtx") jumSun(-25,1,25,1);//for zVtx
    else jumSun(0,1,xmax,1);
    c2->SaveAs(Form("pdf/compareBtwRuns_ratio_refrun26%s_%s%s.pdf",run[0].data(),var.data(),cap));
} 

