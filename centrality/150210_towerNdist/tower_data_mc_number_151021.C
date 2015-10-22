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

void tower_data_mc_number_151021()
{

    const TCut threCut = "abs(eta)>2.87 && abs(eta)<5.2 && et>1.4";
    const int Ncut = 5;
    int col[] = {1,2,4,6,8,28,46,41};
    int marker[] = {20,22,29,33,34};
    string FilterName[Ncut] = {"no evt filter","pprimaryVertexFilter","phltPixelClusterShapeFilter", "phfCoincFilter3", "pcollisionEventSelection"};

    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    // ===================================================================================
    // Get Trees from data & mc files.
    // ===================================================================================

    TFile *f[2]; // 0:data, 1:MC

    TTree *t_evt[2]; 
    TTree *t_skim[2];
    TTree *t_hlt[2];
    TTree *t_recTower[2];
    TTree *t_pfTower[2];
    int Entries[2];   
    for(int i=0;i<2;i++){ 
        if(i==0) f[0] = TFile::Open("root://cluster142.knu.ac.kr//store/user/ygo/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611.root");
        //if(i==0) f[0] = new TFile("root://cluster142.knu.ac.kr//store/user/ygo/PbPb_minbias_data_2760_HIRun2011-14Mar2014-v2_run181611.root");
        else f[1] = TFile::Open("root://cluster142.knu.ac.kr//store/user/ygo/officialMC_Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV.root");
        t_evt[i] = (TTree*) f[i] -> Get("hiEvtAnalyzer/HiTree");
        t_skim[i] = (TTree*) f[i] -> Get("skimanalysis/HltTree");
        t_hlt[i] = (TTree*) f[i] -> Get("hltanalysis/HltTree");
        t_recTower[i] = (TTree*) f[i] -> Get("rechitanalyzer/tower");
        t_pfTower[i] = (TTree*) f[i] -> Get("pfTowers/tower");
        t_recTower[i] -> AddFriend(t_hlt[i]);
        t_recTower[i] -> AddFriend(t_evt[i]);
        t_recTower[i] -> AddFriend(t_skim[i]);
        t_pfTower[i] -> AddFriend(t_hlt[i]);
        t_pfTower[i] -> AddFriend(t_evt[i]);
        t_pfTower[i] -> AddFriend(t_skim[i]);
        //cout << t_evt << ", " << t_skim << ", " << t_hlt << ", " << t_recTower << ", " << t_pfTower << ", " << endl;
        Entries[i] = t_recTower[i] -> GetEntries();
    }

    // ===============================================================================================
    // [et] Define et histograms (data/mc , pf tree/rechit tree) 
    // ===============================================================================================
    cout << "LET'S DEFINE HISTOGRAMS" << endl;

    TH1F* h1F_rec_nTower[3][10]; //[0:data, 1:mc, 2:ratio=data/mc][condition]

    for(int i=0; i<3; i++){
        for(int j=0; j<Ncut; j++){
            h1F_rec_nTower[i][j] = new TH1F(Form("rec_nTower_sample%d_cut%d",i,j), ";# of HF towers;Events",100,0,1000);
            h1F_rec_nTower[i][j] -> SetMarkerStyle(marker[i]);
            h1F_rec_nTower[i][j] -> SetMarkerSize(0.9);
            h1F_rec_nTower[i][j] -> SetMarkerColor(col[i]);
            h1F_rec_nTower[i][j] -> SetLineColor(col[i]);
            h1F_rec_nTower[i][j] -> SetLabelSize(0.03);
            //ratio_rec_n[i] -> SetAxisRange(0.5,1.5,"Y");
        }
    }

    // ===============================================================================================
    // [rechit] Fill histogram by using <Draw> function in TTree.
    // ===============================================================================================
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;
    //********************************************************
    //cut define
    TString towerCut = "abs(tower.eta) > 2.87 && tower.et > 0.5";
    TCut cut="";
    TCut totcut[Ncut];
    totcut[0] = cut;
    totcut[1] = cut&& Form("pprimaryVertexFilter==1");
    //if(isAOD) totcut[2] = cut&& Form("pclusterCompatibilityFilter==1");
    totcut[2] = cut&& Form("phltPixelClusterShapeFilter==1");
    totcut[3] = cut&& Form("phfCoincFilter3==1");
    totcut[4] = cut&& Form("pcollisionEventSelection==1");
    //********************************************************
    for(int i=0; i<2; i++){
        for(int j=0; j<Ncut; j++){
            cout << "i = " << i << ", j = " << j << endl;
            //if(i==0) t_recTower[i]->Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.et > 0.5)>>+%s",h1F_rec_nTower[i][j]->GetName()));
            if(i==0) t_recTower[i]->Draw(Form("Sum$(%s)>>+%s",towerCut.Data(),h1F_rec_nTower[i][j]->GetName()), totcut[j] && "HLT_HIMinBiasHfOrBSC_v1==1");
            else if(i==1) t_recTower[i]->Draw(Form("Sum$(%s)>>+%s",towerCut.Data(),h1F_rec_nTower[i][j]->GetName()), totcut[j]);
            cout << "i = " << i << ", j = " << j << endl;
            h1F_rec_nTower[i][j] = (TH1F*)gDirectory->Get(h1F_rec_nTower[i][j]->GetName());
        }
    }

    // ===============================================================================================
    // [# of HF tower] Normalization!!!!!!! important
    // ===============================================================================================
 
    for(int j=0; j<Ncut; j++){
        int cutBinFrom = h1F_rec_nTower[0][j]->FindBin(200); 
        int cutBinTo = h1F_rec_nTower[0][j]->FindBin(700); 

        h1F_rec_nTower[1][j]->Scale(h1F_rec_nTower[0][j]->Integral(cutBinFrom,cutBinTo)/h1F_rec_nTower[1][j]->Integral(cutBinFrom,cutBinTo));
        h1F_rec_nTower[2][j]->Divide(h1F_rec_nTower[0][j],h1F_rec_nTower[1][j]);

        //data_rec_n[i] -> Scale(mc_rec_n[i]->Integral(cutBinFrom,cutBinTo)/data_rec_n[i]->Integral(cutBinFrom,cutBinTo));
        //ratio_rec_n[i] -> Divide(data_rec_n[i],mc_rec_n[i]);
        //cout << "total integral of MC in the whole range: " << mc_rec_n[i]->Integral()<< endl;
        //cout << "total integral of MC in the cutted range : " << mc_rec_n[i]->Integral(cutBinFrom,cutBinTo)<< endl;
    }


    // ===============================================================================================
    // [# of HF tower] Draw histograms in Canvas. 
    // ===============================================================================================
    cout << "LET'S DRAW HISTOGRAMS IN CANVAS" << endl;
    TCanvas *c_tot = new TCanvas("c_tot","c_tot", 1200,300);
    makeMultiPanelCanvas(c_tot,Ncut,1,0.0,0.0,0.2,0.15,0.02);
    for(int j=0; j<Ncut; j++){
        c_tot->cd(j+1);
        for(int i=0; i<2; i++){
            if(i==0) h1F_rec_nTower[i][j]->DrawCopy();
            else if(i==1) {
                h1F_rec_nTower[i][j]->DrawCopy("hist same");
                drawText(Form("%s",totcut[j].GetTitle()),0.46,0.78);
            }
            if(i==0 && j==0) drawText("rechitTowers", 0.46,0.88); 
        }
        cleverRange(h1F_rec_nTower[0][j],h1F_rec_nTower[1][j]);
    }

    TCanvas *c_filterRate[2];
    for(int i=0;i<2;i++){
        c_filterRate[i] = new TCanvas(Form("c_filterRate_sample%d",i), "c_filterRate", 400,800);
    //    makeMultiPanelCanvas(c_filterRate[i],1,2,0.0,0.0,0.2,0.15,0.02);
        c_filterRate[i]->Divide(1,2);
    }
    for(int i=0; i<2; i++){
        c_filterRate[i]->cd(1);
        gPad->SetLogy();
        for(int j=0; j<Ncut; j++){
            if(j==0) h1F_rec_nTower[i][j]->Draw("ep"); 
            else h1F_rec_nTower[i][j]->Draw("same ep"); 
            h1F_rec_nTower[i][j]->GetXaxis()->SetRangeUser(0.,250.);
        }
    }

    TH1F* h1F_FilterRate[2][Ncut];//[0:data,1:MC][filter]
    for(int i=0; i<2; i++){
        c_filterRate[i]->cd(2);
        gPad->SetLogy();
        for(int j=1; j<Ncut; j++){
            h1F_FilterRate[i][j]= (TH1F*)h1F_rec_nTower[1][0]->Clone(Form("hFilterRate_sample%d_filter%d",i,j));
            h1F_FilterRate[i][j]->GetYaxis()->SetTitle("Filter Efficiency");
            h1F_FilterRate[i][j]->GetXaxis()->SetRangeUser(0.,250.);
            h1F_FilterRate[i][j]->Divide(h1F_rec_nTower[i][j],h1F_rec_nTower[i][0],1,1,"B");
            h1F_FilterRate[i][j]->SetMarkerColor(col[j]);
            h1F_FilterRate[i][j]->SetMarkerStyle(marker[j]);
            if(j==1) h1F_FilterRate[i][j]->DrawCopy();
            else h1F_FilterRate[i][j]->DrawCopy("same");
        }
    }
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    legStyle(leg);
//    leg->SetFillColor(0);
//    leg->SetTextFont(43);
//    leg->SetTextSize(15);
    leg->SetHeader("rechitTowers MC");
    for(int j=1; j<Ncut; j++){
        leg->AddEntry(h1F_rec_nTower[1][j],Form("%s",FilterName[1].c_str()),"p");
    }
    leg->Draw();

    c_filterRate[0] -> SaveAs("png/rectower_filterRate_data_nTower.png");
    c_filterRate[1] -> SaveAs("png/rectower_filterRate_mc_nTower.png");
    c_tot-> SaveAs("png/rectower_data_mc_nTower_tot.png");
    TFile* outFile = new TFile("recTower_n.root","recreate");
    outFile->cd();
    for(int j=0; j<Ncut; j++){
        for(int i=0; i<3; i++){
            h1F_rec_nTower[i][j]->Write();
            if(j!=0&&i!=2) h1F_FilterRate[i][j]->Write();
        }
    }
    outFile->Close();
} 
