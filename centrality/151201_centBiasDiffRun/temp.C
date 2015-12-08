// comparison centrality objects between runs 
// Author : Yeonju Go
// 30 Nov 2015

//private setup
#include "../../yjUtility.h"

const int colhere[] = {1,2,kViolet,4,kOrange,kGreen};
//const int colhere[] = {1,4,9,6,46,kOrange,kViolet,kOrange+10};
//const int colhere[] = {1,2,4,9,6,46,kOrange,kViolet,kOrange+10};
void draw_compare_runs_1D(string var="hiHF");
void temp(){
    TStopwatch timer;
    timer.Start();
    //draw_compare_runs_1D("zVtx");//change many things !! 
    draw_compare_runs_1D("xVtx");//change many things !! 
    //draw_compare_runs_1D("yVtx");//change many things !! 
    //draw_compare_runs_1D("hiHF");
    //draw_compare_runs_1D("hiBin");
    //draw_compare_runs_1D("hiNpix");
    //draw_compare_runs_1D("hiNtracks");
#if 0
    draw_compare_runs_1D("sumPtVtx");
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
    const int Nrun = 6;
    string run[] = {"2816","2811","2695","3233","2620","2656"};
    //string run[] = {"2816","3233"};
    //string run[] = {"2816","2620","2656","2694","2811","2695"};
    //string run[] = {"2695","2548","2620","2656","2694","2811","2816","3031"};
    //string run[] = {"2695","2548","2620","2656","2694","2811","2816","3033"};
    //string run[] = {"695","548","620","656","694","811","816"};
    //const char* run[] = {"695","548","620","656","694"};
    //const char* run[] = {"548","620","656","694","695"};
    //const char* trigCut = "HLT_HIL1MinimumBiasHF1AND_v1";
    TCut centCut = "hiHF>220";

    TCut trigCut = "hiHF>220 && HLT_HIL1MinimumBiasHF1AND_v1 && pcollisionEventSelection";
    TCut trigCutFor3233 = "hiHF>220 && HLT_HIL1MinimumBiasHF2AND_part2_v1 && pcollisionEventSelection";
    //TCut trigCut = "hiBin<150 && HLT_HIL1MinimumBiasHF1AND_v1 && pprimaryVertexFilter && phfCoincFilter3";
    //TCut trigCutFor3233 = "hiBin<150 && HLT_HIL1MinimumBiasHF2AND_part2_v1 && pprimaryVertexFilter && phfCoincFilter3 && lumi>300";
    const char* cap= "compare_2816_3233";

    double trkPtBin[] = {0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0,12.0,16.0,20.0};
    int NtrkPtBin = sizeof(trkPtBin)/sizeof(double)-1;

    ////////////////////// for express
    //const char* dir = "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/";
    //const int Nrun = 8;
    //const char* run[] = {"768","548-v6", "620-v6", "695","703","726","811","837"};
    //const char* run[] = {"548-v6", "620-v6", "656","694","695","697","703","726","768","777","784","817","811","816","818","819","834","836","837"};
    //const char* trigCut = "HLT_HIL1MinimumBiasHF1ANDExpress_v1 && pprimaryVertexFilter && phfCoincFilter3";

    int nBin = 100;
    if(var=="hiBin") nBin=150;
    double xmax = 6000;
    if(var == "hiNpix") xmax = 50000;
    if(var == "hiNtracks" || var == "sumPtVtx") xmax = 3500;
    if(var == "hiBin") xmax = 150;
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
        //if(run[i]=="3233") nEvents[i] = t1->GetEntries(trigCutFor3233);
        //else nEvents[i] = t1->GetEntries(trigCut);
       // cout << run[i].data() << " # of events : " << nEvents[i] << endl;
        if(var == "zVtx") h[i] = new TH1D(Form("h%d",i), Form(";%s (cm);Event Fraction",var.data()), nBin, -25, 25); // for zVtx 
        else if(var == "yVtx"|| var=="xVtx") h[i] = new TH1D(Form("h%d",i), Form(";%s (cm);Event Fraction",var.data()), nBin, 0.06, 0.12); // for zVtx 
        else h[i] = new TH1D(Form("h%d",i), Form(";%s;Event Fraction",var.data()), nBin, 0, xmax);
        if(run[i]=="3233") t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCutFor3233);
        else t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCut);
    }

    TCanvas* c1 = new TCanvas("c2","",500,500);
    if(var != "hiBin") c1->SetLogy();
    for(int i=0; i<Nrun; ++i){
        //h[i]->Scale(1./nEvents[i]); // this is for all!! 
        h[i]->Scale(1./h[i]->Integral("width"));
        SetHistColor(h[i],colhere[i]);
        if(i==0) h[i]->DrawCopy("hist");
        else h[i]->DrawCopy("hist same");
        l1->AddEntry(h[i], Form("Run 26%s",run[i].data()),"pl"); 
    }
    l1->Draw();
    c1->SaveAs(Form("pdf/compareBtwRuns_run26%s_%s%s.pdf",run[0].data(),var.data(),cap));

    h[0]->Rebin(2);
    TCanvas* c2 = new TCanvas("c3","", 500, 500); 
    for(int i=1; i<Nrun; ++i){
        h[i]->Rebin(2);
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

