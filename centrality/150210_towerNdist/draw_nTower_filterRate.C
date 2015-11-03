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

void draw_nTower_filterRate(float etThr=0.5)
{
    const int Ncut = 5;
    int col[] = {1,2,4,6,8,28,46,41};
    int marker[] = {20,22,29,33,34};
    TString FilterName[Ncut] = {"no filter","pprimaryVertexFilter","phltPixelClusterShapeFilter", "phfCoincFilter3", "pcollisionEventSelection"};

    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TString fname = Form("recTower_n_etThr%.1f.root", etThr);
    TFile* fin = new TFile(fname.Data());

    // ===============================================================================================
    // Get nTower histograms
    // ===============================================================================================

    TH1F* h1F_rec_nTower[2][10]; //[0:data, 1:mc, 2:ratio=data/mc][condition]
    TH1F* h1F_FilterRate[2][10]; //[0:data, 1:mc][ratio]

    for(int j=0; j<Ncut; j++){
        for(int i=0; i<2; i++){
            h1F_rec_nTower[i][j] = (TH1F*) fin->Get(Form("rec_nTower_sample%d_filter%d",i,j));
            h1F_rec_nTower[i][j] -> SetMarkerStyle(marker[j]);
            h1F_rec_nTower[i][j] -> SetMarkerSize(0.9);
            h1F_rec_nTower[i][j] -> SetMarkerColor(col[j]);
            h1F_rec_nTower[i][j] -> SetLineColor(col[j]);
            h1F_rec_nTower[i][j] -> SetLabelSize(0.03);
            h1F_FilterRate[i][j] = (TH1F*) h1F_rec_nTower[i][j]->Clone(Form("rec_nTower_sample%d_filterRate%d",i,j));
            h1F_FilterRate[i][j]->SetName(";# of HF towers;Filter Rate");
            //if(j!=0) h1F_FilterRate[i][j]->Divide(h1F_rec_nTower[i][j],h1F_rec_nTower[i][0]);
        }
    }

    TLegend *leg = new TLegend(0.2,0.55,0.7,0.8);
    legStyle(leg);
    //leg->SetHeader("rechitTowers");
    for(int j=0; j<Ncut; j++){
        leg->AddEntry(h1F_rec_nTower[1][j],Form("%s",FilterName[j].Data()),"p");
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
            //h1F_rec_nTower[i][j]->GetXaxis()->SetRangeUser(0.,250.);
        }
        if(i==0) drawText(Form("DATA, tower et>%.1f",etThr),0.46,0.80);
        else if(i==1) drawText(Form("MC, tower et>%.1f",etThr),0.46,0.80);
        leg->Draw("same");
    }

    TLegend *legfR[2];
    for(int i=0; i<2; i++){
        legfR[i] = new TLegend(0.2,0.25,0.9,0.55);
        legStyle(legfR[i]);
    }
    float dy = 0.06;
    for(int i=0; i<2; i++){
        c_filterRate[i]->cd(2);
        gPad->SetLogx();
        for(int j=1; j<Ncut; j++){
            h1F_FilterRate[i][j]->GetYaxis()->SetTitle("Filter Rate");
            h1F_FilterRate[i][j]->Divide(h1F_rec_nTower[i][j],h1F_rec_nTower[i][0],1,1,"B");
            h1F_FilterRate[i][j]->SetMarkerColor(col[j]);
            h1F_FilterRate[i][j]->SetMarkerStyle(marker[j]);
            if(etThr==0.5) h1F_FilterRate[i][j]->GetYaxis()->SetRangeUser(0.9,1.01);
            else if(etThr==1.0) h1F_FilterRate[i][j]->GetYaxis()->SetRangeUser(0.994,1.001);
            else if(etThr==1.5) h1F_FilterRate[i][j]->GetYaxis()->SetRangeUser(0.994,1.001);
            if(j==1) {h1F_FilterRate[i][j]->DrawCopy();}
            else h1F_FilterRate[i][j]->DrawCopy("same");
        }
        jumSun(0,1,1000,1);
        float inte_nofilter = h1F_rec_nTower[i][0]->Integral();
        for(int j=1; j<Ncut; j++){
            float inte_filter = h1F_rec_nTower[i][j]->Integral();
            legfR[i]->AddEntry(h1F_rec_nTower[1][j],Form("%s Rate = %.6f",FilterName[j].Data(),inte_filter/inte_nofilter),"p");
            //drawText(Form("%s Rate = %.4f",FilterName[j].Data(), inte_filter/inte_nofilter),0.3,0.6-dy*j);
        }
        c_filterRate[i]->cd(2);
        legfR[i]->Draw("same");
    }

    c_filterRate[0] -> SaveAs(Form("png/nTower_filterRate_data_etThr%.1f.png",etThr));
    c_filterRate[1] -> SaveAs(Form("png/nTower_filterRate_mc_etThr%.1f.png",etThr));

    c_filterRate[0] -> SaveAs(Form("pdf/nTower_filterRate_data_etThr%.1f.pdf",etThr));
    c_filterRate[1] -> SaveAs(Form("pdf/nTower_filterRate_mc_etThr%.1f.pdf",etThr));
} 

