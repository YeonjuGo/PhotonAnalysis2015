#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLatex.h"
#include "TCut.h"
#include "stdio.h"
#include "../yjUtility.h"

const double Pi = TMath::Pi();
const int col[] = {1,2,4};
void drawSame(TH1* h1=0, TH1* h2=0, TH1* h3=0, const char* mode="AllQCDPhoton30", bool isLogY=1);
void temp(const char* mode= "AllQCDPhoton30")
{
    TH1::SetDefaultSumw2();
    gStyle->SetLabelSize(0.03,"Y");
    gStyle -> SetTitleYSize(0.05);
    // gStyle -> SetTitleXSize(0.06);
    const int nfile = 3;
    TFile *f[nfile];
    TTree* t[nfile];
    f[0] = new TFile(Form("files/step3_%s_multiFit25ns.root",mode));
    t[0] = (TTree*) f[0]->Get("Events");
    f[1] = new TFile(Form("files/step3_%s_53XecalLocalReco.root",mode));
    t[1] = (TTree*) f[1]->Get("Events");
    f[2] = new TFile(Form("files/step3_%s_multifit50ns.root",mode));
    t[2] = (TTree*) f[2]->Get("Events");

    TH1D* e_EB[nfile];
    TH1D* time_EB[nfile];
    TH1D* e_EE[nfile];
    TH1D* time_EE[nfile];

    TH1D* phoEt_EE[nfile];
    TH1D* phoEt_EB[nfile];
    TH1D* phoPhi_EE[nfile];
    TH1D* phoPhi_EB[nfile];

    for(int i=0; i<nfile; i++){
        e_EB[i] = new TH1D(Form("rechit_energy_EB_sample%d",i),";EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO energy;",50,0,70);
        time_EB[i] = new TH1D(Form("rechit_time_EB_sample%d",i),";EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO time;",50,-150,150);
        e_EE[i] = (TH1D*)e_EB[i]->Clone(Form("rechit_energy_EE_sample%d",i));
        e_EE[i]->SetTitle(";EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO energy;");
        time_EE[i] = (TH1D*) time_EB[i]->Clone(Form("rechit_time_EE_sample%d",i));
        time_EE[i]->SetTitle(";EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO time;");

        phoEt_EB[i] = new TH1D(Form("phoEt_EB_sample%d",i),";recoPhotons_photons__RECO barrel phoEt;",110,0,110);
        phoEt_EE[i] = (TH1D*) phoEt_EB[i]->Clone(Form("phoEt_EE_sample%d",i));
        phoEt_EE[i]->SetTitle(";recoPhotons_photons__RECO endcap phoEt;");
        phoPhi_EB[i] = new TH1D(Form("phoPhi_EB_sample%d",i),";recoPhotons_photons__RECO barrel phoPhi;",50,Pi,-Pi);
        phoPhi_EE[i] = (TH1D*) phoPhi_EB[i]->Clone(Form("phoPhi_EE_sample%d",i));
        phoPhi_EE[i]->SetTitle(";recoPhotons_photons__RECO endcap phoPhi;");
    }

    for(int i=0; i<nfile; i++){
        t[i]->Draw(Form("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.obj.obj.energy_>>+%s",e_EB[i]->GetName()));
        t[i]->Draw(Form("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.obj.obj.time_>>%s",time_EB[i]->GetName()));
        t[i]->Draw(Form("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_RECO.obj.obj.energy_>>%s",e_EE[i]->GetName()));
        t[i]->Draw(Form("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_RECO.obj.obj.time_>>%s",time_EE[i]->GetName()));
        t[i]->Draw(Form("recoPhotons_photons__RECO.obj.m_state.p4Polar_.Et()>>%s",phoEt_EB[i]->GetName()));
        t[i]->Draw(Form("recoPhotons_photons__RECO.obj.m_state.p4Polar_.Phi()>>%s",phoPhi_EB[i]->GetName()));
    }

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    // Draw
    drawSame(e_EB[0], e_EB[1], e_EB[2], mode,1);
    drawSame(time_EB[0], time_EB[1], time_EB[2], mode,0);
    drawSame(e_EE[0], e_EE[1], e_EE[2], mode,1);
    drawSame(time_EE[0], time_EE[1], time_EE[2], mode,0);
    drawSame(phoEt_EB[0], phoEt_EB[1], phoEt_EB[2], mode,0);
    drawSame(phoPhi_EB[0], phoPhi_EB[1], phoPhi_EB[2], mode,0);
}

void drawSame(TH1* h1, TH1* h2, TH1* h3, const char* mode, bool isLogY){

    TCanvas* c=  new TCanvas(Form("c_%s",h1->GetName()),"", 400,400);
    h1->Sumw2();
    h2->Sumw2();
    h3->Sumw2();
    //  h1->Scale( 1. / h1->Integral());
    //  h2->Scale( 1. / h2->Integral());

    // default plotting options
    h1->SetMarkerStyle(kFullSquare);
    h2->SetMarkerStyle(kFullCircle);
    h3->SetMarkerStyle(kFullCircle);

    h1->SetMarkerColor(col[0]);
    h2->SetMarkerColor(col[1]);
    h3->SetMarkerColor(col[2]);

    h3->SetMarkerSize(h3->GetMarkerSize()*0.8);     // to distinguish between points when they overlap

    // set Y-axis ranges
    // default Y-axis range is that of h1
    // make sure that both plots will not run out of y-axis

    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
    legStyle(leg);
    leg->AddEntry(h1, "75X");
    leg->AddEntry(h2, "75X with 53X ecalReco");
    leg->AddEntry(h3, "75X with 50 ns");
    /*
       h1->DrawCopy("p");
       h2->DrawCopy("p same");
       h3->DrawCopy("p same");
     */
    h1->Draw("p");
    //    leg->Draw("same");
    h2->Draw("same");
    h3->Draw("same");
    if(isLogY) gPad->SetLogy();
    c -> SaveAs(Form("pdf/%s_%s.pdf",mode,h1->GetName())); 
}

int main(){
    temp("AllQCDPhoton30");
    return 0;
}
