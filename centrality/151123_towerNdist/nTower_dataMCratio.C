// Author Yeonju Go
// last modification : 2015/01/27 
// 
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
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
#include "TLatex.h"
#include "stdio.h"
#include "../../yjUtility.h"

float normHist(TH1* hNom=0, TH1* hDen=0, TH1* hRatio=0, double cut_i=700, double cut_f=900);
double makeHist_nTower(float etThr=0.5, float eThr=0, float norm_i=300, float norm_f=400);

void nTower_dataMCratio(){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);

    double thr[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
    int nbin = sizeof(thr)/sizeof(double);
    double thrArr[nbin+1];
    for(int i=0;i<nbin+1;i++){
        if(i!=nbin) thrArr[i] = thr[i]-0.25;
        else thrArr[i] = thr[i]+0.25;
    }
    
    TH1D* h1D_e = new TH1D("h1D_e_ratio",";(# of HF tower) Energy Threshold (GeV);DATA/MC Ratio", nbin,thrArr);
    for(int i=0;i<nbin;i++){
        cout << "START energy threshold : " << thr[i] << endl;
        double cont = makeHist_nTower(0,thr[i],300,400); 
        h1D_e->SetBinContent(i,cont);
    }

    TCanvas* c_e = new TCanvas();
    h1D_e->Draw("hist");
    c_e->SaveAs("png/eThr_vs_ratio_distribution.png");

}


float normHist(TH1* hNom, TH1* hDen, TH1* hRatio, double cut_i, double cut_f){
    int cutBinFrom = hDen->FindBin(cut_i);
    int cutBinTo = hDen->FindBin(cut_f);

    TH1::SetDefaultSumw2();
    hDen->Scale(hNom->Integral(cutBinFrom,cutBinTo)/hDen->Integral(cutBinFrom,cutBinTo));
    hRatio->Divide(hNom,hDen);

    //cout << "total integral of DATA in the whole range: " << hData->Integral()<< endl;
    //cout << "total integral of MC in the whole range: " << hMC->Integral()<< endl;
    //cout << "full efficiency in the whole range: " <<  hData->Integral()/hMC->Integral()<< endl;

    return hNom->Integral()/hDen->Integral();
}

double makeHist_nTower(float etThr, float eThr, float norm_i, float norm_f)
{
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);
    TH1::SetDefaultSumw2();

    TCut dataCut("HLT_HIL1MinimumBiasHF1AND_v1 && phfCoincFilter3");
    TCut mcCut("");

    // ===================================================================================
    // Get Trees from data & mc files.
    // ===================================================================================
    TFile *f[2]; // 0:data, 1:MC

    TTree *t_skim[2];
    TTree *t_hlt[2];
    TTree *t[2];
    int Entries[2];   
    const char* fname_data="root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HiForest_Streamer_run262315.root";
    const char* fname_mc="/afs/cern.ch/work/y/ygo/public/centrality/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_merged_forest_6.root";

    for(int i=0;i<2;i++){ 
        if(i==0) f[0] = TFile::Open(fname_data);
        else f[1] = TFile::Open(fname_mc);
        t_skim[i] = (TTree*) f[i] -> Get("skimanalysis/HltTree");
        t_hlt[i] = (TTree*) f[i] -> Get("hltanalysis/HltTree");
        t[i] = (TTree*) f[i] -> Get("rechitanalyzer/tower");
        t[i] -> AddFriend(t_hlt[i]);
        t[i] -> AddFriend(t_skim[i]);
        if(i==0) Entries[i] = t[i] -> GetEntries(dataCut);
        else Entries[i] = t[i] -> GetEntries();
    }

    // ===============================================================================================
    // [et] Define et histograms (data/mc , pf tree/rechit tree) 
    // ===============================================================================================
    cout << "LET'S DEFINE HISTOGRAMS" << endl;
    TH1F* h1F[2]; //[0:data, 1:mc]
    for(int i=0; i<2; i++){
        h1F[i] = new TH1F(Form("h1F_sample%d",i), ";# of HF towers;Events",200,0,1000);
        h1F[i] -> SetMarkerStyle(marker[i]);
        h1F[i] -> SetMarkerSize(0.9);
        if(i==0) SetHistColor(h1F[i],2);
        else SetHistColor(h1F[i],1);
    }

    // ===============================================================================================
    // [rechit] Fill histogram by using <Draw> function in TTree.
    // ===============================================================================================
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;
    //********************************************************
    const char* towerCut = Form("abs(tower.eta) > 2.87 && abs(tower.eta)<5.205 && tower.et*cosh(eta) > %.1f && tower.et > %.1f",eThr, etThr);
    
    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();

    for(int i=0; i<2; i++){
        std::clock_t    startT_draw, endT_draw;
        startT_draw = std::clock();
        if(i==0) t[i]->Draw(Form("Sum$(%s)>>+%s",towerCut,h1F[i]->GetName()), dataCut);
        else if(i==1) t[i]->Draw(Form("Sum$(%s)>>+%s",towerCut,h1F[i]->GetName()), mcCut);
        h1F[i] = (TH1F*)gDirectory->Get(h1F[i]->GetName());
        endT_draw = std::clock();
        std::cout.precision(6);      // get back to default precision
        std::cout << "DRAW finished in             : " << (endT_draw - startT_draw) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    }

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited loop" << std::endl;

    TH1F* h1F_ratio;
    h1F_ratio = (TH1F*) h1F[1]->Clone(Form("h1F_ratio"));
    h1F_ratio->SetName(";# of HF towers;DATA/MC");

    // ===============================================================================================
    // [# of HF tower] Normalization!!!!!!! important
    // ===============================================================================================
    float ratio = normHist(h1F[0],h1F[1], h1F_ratio, norm_i, norm_f);

    // ===============================================================================================
    // [# of HF tower] Draw histograms in Canvas.
    // ===============================================================================================
    cout << "LET'S DRAW HISTOGRAMS IN CANVAS" << endl;
    TCanvas *c = new TCanvas("c","c", 500,1000);
    c->Divide(2,1);
    c->cd(1);
    gPad->SetLogy();
    h1F[0]->DrawCopy();
    h1F[1]->DrawCopy("hist same");
    drawText(Form("tower et > %.1f",etThr), 0.46,0.80);
    drawText(Form("tower e > %.1f",eThr), 0.46,0.88);
    Double_t range = cleverRange(h1F[0],h1F[1]);
    jumSun(norm_i, 0.000001, norm_i, range);
    jumSun(norm_f, 0.000001, norm_f, range);
    c->cd(2);    
    h1F_ratio->GetYaxis()->SetTitle("DATA/MC");
    h1F_ratio->DrawCopy();
    drawText(Form("DATA/MC = %.3f", ratio), 0.26,0.88);
    c->SaveAs(Form("png/nTower_eThr%.1f_etThr%.1f_normRange%dto%d.png",eThr,etThr,(int)norm_i,(int)norm_f));

    TFile* outFile = new TFile(Form("recTower_n_eThr%.1f_etThr%.1f.root",eThr,etThr),"recreate");
    outFile->cd();
    h1F[0]->Write();
    h1F[1]->Write();
    h1F_ratio->Write();
    outFile->Close();
/*
    delete h1F_ratio; 
    delete h1F[0]; 
    delete h1F[1]; 
    delete f[0];
    delete f[1];
    delete t_skim[0];
    delete t_skim[1];
    delete t_hlt[0];
    delete t_hlt[1];
    delete t[0];
    delete t[1];
  */  
    return ratio;
} 
