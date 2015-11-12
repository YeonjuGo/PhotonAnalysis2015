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

    const char* infile =  Form("histfiles/hist_noiseCompDiffTrigger.root");
    const char* trig[]={"HLT_HIZeroBias_v1", "HLT_HIMinBiasZDC_Calo_v1", "HLT_HIMinBiasHfOrBSC_v1", "!HLT_HIMinBiasHfOrBSC_v1 && HLT_HIZeroBias_v1"};
    const int nTrig = 4;

    TFile * inf = new TFile(infile);

    TH1D* e[nTrig];
    TH1D* et[nTrig];
    TH1D* e_tower[nTrig];
    TH1D* et_tower[nTrig];
    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
    for(int i = 0; i < nTrig; ++i){
        e[i] = (TH1D*) inf->Get(Form("e%d",i));
        et[i] = (TH1D*) inf->Get(Form("et%d",i));
        e_tower[i] = (TH1D*) inf->Get(Form("e_tower%d",i));
        et_tower[i] = (TH1D*) inf->Get(Form("et_tower%d",i));
        if(i==0) lge->AddEntry(e[i], trig[i]);
        e[i]->SetMarkerColor(col[i]); 
        et[i]->SetMarkerColor(col[i]);
        e_tower[i]->SetMarkerColor(col[i]);
        et_tower[i]->SetMarkerColor(col[i]);
    }

    TCanvas* c1 = new TCanvas("c1","",400,400);
    for(int i = 0; i < nTrig; ++i){
        if(i==0) e[i]->Draw();
        else e[i]->Draw("same");
    }
    leg->Draw("same");
    c1->SaveAs("pdf/hfe_noiseCompDiffTrigger.pdf");
}

