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

float mean(float data[], int n);
float standard_deviation(float data[], int n);
float normHist(TH1* hMC=0, TH1* hData=0, TH1* hRatio=0, double cut_i=700, double cut_f=900){
    int cutBinFrom = hMC->FindBin(cut_i); 
    int cutBinTo = hMC->FindBin(cut_f); 

    hMC->Scale(hData->Integral(cutBinFrom,cutBinTo)/hMC->Integral(cutBinFrom,cutBinTo));
    hRatio->Divide(hData,hMC);

    //cout << "total integral of DATA in the whole range: " << hData->Integral()<< endl;
    //cout << "total integral of MC in the whole range: " << hMC->Integral()<< endl;
    //cout << "full efficiency in the whole range: " <<  hData->Integral()/hMC->Integral()<< endl; 

    return hData->Integral()/hMC->Integral();
}

void draw_nTower(float etThr=0.0, float eThr=5.0, int N=50, float norm_i=400, float norm_f=1200)
{
    const int Ncut = 5;
    int col[] = {1,2,4,6,8,28,46,41};
    int marker[] = {20,22,29,33,34};
    TString FilterName[Ncut] = {"no filter","pprimaryVertexFilter","phltPixelClusterShapeFilter", "phfCoincFilter3", "pcollisionEventSelection"};

    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TString fname = Form("histfiles_5320/sumHF_etThr%.1f_eThr%.1f.root", etThr,eThr);
    TFile* fin = new TFile(fname.Data());

    // ===============================================================================================
    // [et] Define et histograms (data/mc , pf tree/rechit tree) 
    // ===============================================================================================

    TH1F* h1F_sumHF_data[10]; //[condition]
    TH1F* h1F_sumHF_mc; 
    TH1F* h1F_sumHF_ratio[10]; //[condition]

    for(int j=0; j<Ncut; j++){
        h1F_sumHF_data[j] = (TH1F*) fin->Get(Form("sumHF_data_filter%d",j));
        h1F_sumHF_data[j] -> Rebin(5); 
        h1F_sumHF_data[j] -> SetMarkerStyle(marker[j]);
        h1F_sumHF_data[j] -> SetMarkerSize(0.9);
        h1F_sumHF_data[j] -> SetMarkerColor(col[j]);
        h1F_sumHF_data[j] -> SetLineColor(col[j]);
        h1F_sumHF_data[j] -> SetLabelSize(0.03);
        h1F_sumHF_ratio[j] = (TH1F*) h1F_sumHF_data[j]->Clone(Form("sumHF_ratio_filter%d",j));
        h1F_sumHF_ratio[j]->SetName(";# of HF towers;DATA/MC");
    }
    h1F_sumHF_mc = (TH1F*) fin->Get(Form("sumHF_mc"));
    h1F_sumHF_mc -> Rebin(5); 
    // ===============================================================================================
    // [# of HF tower] Normalization!!!!!!! important
    // ===============================================================================================
    float effi_mean[10];
    float effi_stdv[10];
    for(int j=0; j<Ncut; j++){
        float effi[N];
        for(int jn=0; jn<N; jn++){
            effi[jn] = normHist(h1F_sumHF_mc, h1F_sumHF_data[j],h1F_sumHF_ratio[j], norm_i-(N/2.)+jn, norm_f);
        }
        effi_mean[j] = mean(effi, N);
        effi_stdv[j] = standard_deviation(effi, N);
        //cout << "filter " << j << " ::: efficiency mean = " << effi_mean[j] << ", stdv = " << effi_stdv[j] << endl;
    }

    // ===============================================================================================
    // [# of HF tower] Draw histograms in Canvas. 
    // ===============================================================================================
    cout << "LET'S DRAW HISTOGRAMS IN CANVAS" << endl;
    TCanvas *c_tot = new TCanvas("c_tot","c_tot", 1200,600);
    //makeMultiPanelCanvas(c_tot,Ncut,2,0.0,0.0,0.2,0.15,0.02);
    c_tot->Divide(Ncut,2);
    for(int j=0; j<Ncut; j++){
        c_tot->cd(j+1);
        h1F_sumHF_data[j]->Draw();
        h1F_sumHF_mc->Draw("hist same");
        drawText(Form("%s",FilterName[j].Data()),0.30,0.78);
        //drawText(Form("%s",totcut[j].GetTitle()),0.46,0.78); //TCut.GetTitle() works
        if(j==0) drawText("rechitTowers", 0.46,0.88); 
        else if(j==1) {
            if(eThr!=0.0) drawText(Form("tower e > %.1f",etThr), 0.46,0.68); 
            if(etThr!=0.0) drawText(Form("tower et > %.1f",eThr), 0.46,0.60); 
        }
   //     gPad->SetLogy();
        Double_t range = cleverRange(h1F_sumHF_data[j],h1F_sumHF_data[j]);
        onSun(norm_i, 0.000001, norm_i, range);
        onSun(norm_f, 0.000001, norm_f, range);
        jumSun(norm_i-N/2., 0.000001, norm_i-N/2., range);
        jumSun(norm_i+N/2., 0.000001, norm_i+N/2., range);
        c_tot->cd(j+6);
        h1F_sumHF_ratio[j]->GetYaxis()->SetTitle("DATA/MC");
        h1F_sumHF_ratio[j]->DrawCopy();
        drawText(Form("DATA/MC = %.3f +- %.3f", effi_mean[j], effi_stdv[j]), 0.26,0.88);
        cleverRange(h1F_sumHF_data[j],h1F_sumHF_mc);
    }
    c_tot->SaveAs(Form("pdf/nHF_etThr%.1f_eThr%.1f_normRange%dto%d_stdvVar%d.pdf",etThr,eThr,(int)norm_i,(int)norm_f,N));
    
    TCanvas *c33 = new TCanvas("c33","c33", 600,600);
    c33->Divide(2,2);
    //makeMultiPanelCanvas(c33,2,2,0.0,0.0,0.2,0.15,0.02);
    c33->cd(1);
    h1F_sumHF_data[0]->Draw();
    h1F_sumHF_mc->Draw("hist same");
    cleverRange(h1F_sumHF_data[0],h1F_sumHF_mc);
    c33->cd(2);
    h1F_sumHF_mc->Draw("hist");
    h1F_sumHF_data[Ncut-1]->Draw("same");
    cleverRange(h1F_sumHF_data[Ncut-1],h1F_sumHF_mc);
    c33->cd(3);
    h1F_sumHF_ratio[0]->GetYaxis()->SetTitle("DATA/MC");
    h1F_sumHF_ratio[0]->GetYaxis()->SetRangeUser(0,5);
    h1F_sumHF_ratio[0]->DrawCopy();
    jumSun(0,1,2000,1,2);
    drawText(Form("DATA/MC = %.3f", effi_mean[0]), 0.26,0.70);
    //drawText(Form("DATA/MC = %.3f +- %.3f", effi_mean[0], effi_stdv[0]), 0.26,0.88);
    c33->cd(4);
    h1F_sumHF_ratio[Ncut-1]->GetYaxis()->SetTitle("DATA/MC");
    h1F_sumHF_ratio[Ncut-1]->GetYaxis()->SetRangeUser(0,5);
    h1F_sumHF_ratio[Ncut-1]->DrawCopy();
    jumSun(0,1,2000,1,2);
    
    drawText(Form("DATA/MC = %.3f", effi_mean[Ncut-1]), 0.26,0.70);
    //drawText(Form("DATA/MC = %.3f +- %.3f", effi_mean[Ncut-1], effi_stdv[Ncut-1]), 0.26,0.88);

    c33->SaveAs(Form("pdf/nHF_onlyCollEvtFilter_etThr%.1f_eThr%.1f_normRange%dto%d_stdvVar%d.pdf",etThr,eThr,(int)norm_i,(int)norm_f,N));

    TCanvas *c4;
    //= new TCanvas("c4","c4", 600,600);
    c4 = (TCanvas*)c33->DrawClone();
    for(int i=0;i<Ncut;i++){
        c4->cd(i+1);
        gPad->SetLogx();
    }
    c4->SaveAs(Form("pdf/nHF_onlyCollEvtFilter_etThr%.1f_eThr%.1f_normRange%dto%d_stdvVar%d_logx.pdf",etThr,eThr,(int)norm_i,(int)norm_f,N));

}

float mean(float data[], int n)
{
    float mean=0.0;
    int i;
    for(i=0; i<n;++i)
    {
        mean+=data[i];
    }
    mean=mean/n;
    return mean;           
}
float standard_deviation(float data[], int n)
{
    float mean=0.0, sum_deviation=0.0;
    int i;
    for(i=0; i<n;++i)
    {
        mean+=data[i];
    }
    mean=mean/n;
    for(i=0; i<n;++i)
        sum_deviation+=(data[i]-mean)*(data[i]-mean);
    return sqrt(sum_deviation/n);           
}


int main(){
    draw_nTower(0.0, 1.0);
    draw_nTower(0.0, 3.0);
    draw_nTower(0.0, 5.0);
    return 0;
}
