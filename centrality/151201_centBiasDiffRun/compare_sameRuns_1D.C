// Author : Yeonju Go

#include "../../yjUtility.h"

const int colhere[] = {1,2,4,9,6,46,kOrange,kViolet,kOrange+10};
void draw_compare_sameRuns_1D(const char* fprom="", const char* fexp="", string var="", TCut cut="", string run="");
void compare_sameRuns_1D(){
    TStopwatch timer;
    timer.Start();

    const int Nrun = 7;
    string run[] = {"2620","2656","2694","2695","2697","2811","2816"};
    TCut lumiCut[] = {"lumi>=100 && lumi<=368","lumi>=100 && lumi<=220","lumi<=178","lumi>=71 && lumi<=229","lumi<=227","lumi<=63","lumi>=7 && lumi<=447","lumi<=449"};
    string fexp = "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run26";
    string fprom = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIForestExpress_run26"; //for 620,656,694,695,697
    string fprom2 = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIMinimumBias2_run26"; //for 811,816,3035
    for(int i=0; i<Nrun; i++){
        const char* fin_exp = Form("%s%s.root",fexp.data(),run[i].data());
        const char* fin_prom = Form("%s%s.root",fprom.data(),run[i].data());
        if(run[i] == "2837") fin_exp = Form("%s%s_preHBHE.root",fexp.data(),run[i].data());
        if(run[i] == "2811" || run[i] == "2816") fin_prom = Form("%s%s.root",fprom2.data(),run[i].data());
        if(run[i] == "2548" || run[i] == "2620") fin_exp = Form("%s%s-v6.root",fexp.data(),run[i].data());

        string var[] = {"hiBin","hiHF","hiNpix"};
        //string var[] = {"hiBin","hiHF","hiNpix","hiNtracks"};
        int nVar = 3; 
        for(int j=0;j<nVar;j++){
            draw_compare_sameRuns_1D(fin_prom, fin_exp, var[j], lumiCut[i], run[i]);
        }
    }


    timer.Stop();
    cout<<"Macro finished: "<<endl;
    cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
    cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
}

void draw_compare_sameRuns_1D(const char* fprom, const char* fexp, string var, TCut cut, string run){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2(); 
    const char* cap= "";
    int nBin = 50;
    if(var=="hiBin") nBin=50;
    double xmax = 6000;
    if(var == "hiNpix") xmax = 50000;
    if(var == "hiNtracks" || var == "sumPtVtx") xmax = 3500;
    if(var == "hiBin") xmax = 200;
    if(var == "hiEE") xmax = 3000;
    if(var == "hiEB") xmax = 4000;

    int nEvents[2];
    TH1D* h[2];
    TLegend* l1 = new TLegend(0.7,0.65,0.95,0.95,Form("Run 26%s",run.data()));
    legStyle(l1);
    TCut trigCut[2];
    trigCut[0] = cut && "HLT_HIL1MinimumBiasHF1AND_v1 && pprimaryVertexFilter && phfCoincFilter3";
    trigCut[1] = cut && "HLT_HIL1MinimumBiasHF1ANDExpress_v1 && pprimaryVertexFilter && phfCoincFilter3";
    for(int i=0; i<2; ++i){
        TFile* f1;
        if(i==0) f1 = TFile::Open(fprom);
        else f1 = TFile::Open(fexp);
        cout << "Open file : " << f1->GetName() << endl;
        TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
        TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
        TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
        TTree* t1_track= (TTree*) f1 -> Get("anaTrack/trackTree");
        t1->AddFriend(t1_skim);
        t1->AddFriend(t1_hlt);
        t1->AddFriend(t1_track);
        nEvents[i] = t1->GetEntries(trigCut[i]);
        cout << i << " # of events : " << nEvents[i] << endl;
        if(var == "zVtx") h[i] = new TH1D(Form("h%d",i), Form(";%s (cm);Event Fraction",var.data()), nBin, -25, 25); // for zVtx 
        else h[i] = new TH1D(Form("h%d",i), Form(";%s;Event Fraction",var.data()), nBin, 0, xmax);
        t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCut[i]);
    }

    TCanvas* c1 = new TCanvas("c2","",500,500);
    if(var != "hiBin") c1->SetLogy();
    for(int i=0; i<2; ++i){
        h[i]->Scale(1./nEvents[i]); // this is for all!! 
        SetHistColor(h[i],colhere[i]);
        if(i==0) h[i]->DrawCopy("hist");
        else h[i]->DrawCopy("hist same");
    }
    l1->AddEntry(h[0], Form("Private Reco"),"pl"); 
    l1->AddEntry(h[1], Form("Express"),"pl"); 
    l1->Draw();
    c1->SaveAs(Form("pdf/compareSameRuns_run26%s_%s%s.pdf",run.data(),var.data(),cap));

    h[0]->Rebin(2);
    TCanvas* c2 = new TCanvas("c3","", 500, 500); 
    for(int i=1; i<2; ++i){
        h[i]->Rebin(2);
        h[i]->Divide(h[0]);
        h[i]->SetMarkerStyle(20);
        h[i]->SetMarkerSize(0.85);
        h[i]->GetYaxis()->SetRangeUser(0.6,1.3);
        h[i]->SetTitle(Form(";%s;Express / Private Reco",var.data()));
        if(i==1) h[i]->DrawCopy("pe");
        else h[i]->DrawCopy("pe same");
    }
    if(var == "zVtx") jumSun(-25,1,25,1);//for zVtx
    else jumSun(0,1,xmax,1);
    c2->SaveAs(Form("pdf/compareSameRuns_ratio_run26%s_%s%s.pdf",run.data(),var.data(),cap));
} 

