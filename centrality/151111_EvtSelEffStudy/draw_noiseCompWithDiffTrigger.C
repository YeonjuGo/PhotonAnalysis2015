// original macro by Yetkin?
// modified by yeonju
// to get number of the towers?

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
#include "../../yjUtility.h"

void draw_noiseCompWithDiffTrigger(){
    gStyle->SetOptStat(0);
    const char* infile =  Form("histfiles/hist_noiseCompDiffTrigger.root");
    const char* trig[]={"HLT_HIZeroBias_v1", "HLT_HIMinBiasZDC_Calo_v1", "HLT_HIMinBiasHfOrBSC_v1", "!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1"};
    const int nTrig = 4;

    TFile * inf = new TFile(infile);

    TH1D* e[nTrig];
    TH1D* et[nTrig];
    TH1D* e_tower[nTrig];
    TH1D* et_tower[nTrig];
    TLegend* leg = new TLegend(0.2,0.7,0.9,0.9);
    legStyle(leg);
    for(int i = 0; i < nTrig; ++i){
        e[i] = (TH1D*) inf->Get(Form("e%d",i));
        et[i] = (TH1D*) inf->Get(Form("et%d",i));
        e_tower[i] = (TH1D*) inf->Get(Form("e_tower%d",i));
        et_tower[i] = (TH1D*) inf->Get(Form("et_tower%d",i));
        leg->AddEntry(e[i], trig[i]);
        SetHistColor(e[i],col[i]); 
        SetHistColor(et[i],col[i]);
        SetHistColor(e_tower[i],col[i]);
        SetHistColor(et_tower[i],col[i]);
   }

    TCanvas* c1 = new TCanvas("c1","",400,400);
    for(int i = 0; i < nTrig; ++i){
        e[i]->GetXaxis()->SetRangeUser(0,70);
        e[i]->GetYaxis()->SetRangeUser(0.1,e[i]->GetMaximum()*10); 
        if(i==0) e[i]->DrawCopy();
        else e[i]->DrawCopy("same");
    }
    gPad->SetLogy(); 
    leg->Draw();
    c1->SaveAs("pdf/hfe_noiseCompDiffTrigger.pdf"); 
    TCanvas* c2 = new TCanvas("c2","",400,400);
    for(int i = 0; i < nTrig; ++i){
        if(i==0) {e[i]->GetXaxis()->SetRangeUser(0,15);e[i]->DrawCopy();}
        else e[i]->DrawCopy("same");
    }
    gPad->SetLogy(); 
    leg->Draw();
    c2->SaveAs("pdf/hfe_noiseCompDiffTrigger_zoom.pdf");
/* 
    TCanvas* c3 = new TCanvas("c3","",400,400);
    for(int i = 0; i < nTrig; ++i){
        et[i]->GetXaxis()->SetRangeUser(0,10);
        et[i]->GetYaxis()->SetRangeUser(0.1,et[i]->GetMaximum()*10); 
        if(i==0) et[i]->DrawCopy();
        else et[i]->DrawCopy("same");
    }
    gPad->SetLogy(); 
    leg->Draw();
    c1->SaveAs("pdf/hfet_noiseCompDiffTrigger.pdf");   
*/
}

