#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TLine.h>

void plotPhoton(char *infname = "photonTree.root")
{
    TFile *inf = new TFile(infname);
    TTree *t = (TTree*)inf->Get("t");
    
    TCanvas *c= new TCanvas("c","",600,600);
    
    TH2D *h = new TH2D("h","",12,-2.4,2.4,100,0,2);
    
    t->Draw("phoEt/matchedEt:phoEta>>h","matchedEt>20&&phoHoverE<0.1&&matchedCalIsoDR04<5&&hiBin<200");
    
    h->FitSlicesY(0,0,-1,0,"LL m");
    TH1D *hMean = (TH1D*) gDirectory->Get("h_1");
    hMean->Draw();
    hMean->SetTitle("HLT MultiFit 25ns 0-100%");
    hMean->SetXTitle("Photon #eta");
    hMean->SetYTitle("Reco E_{T} / Gen E_{T}");
    hMean->SetAxisRange(0.9,1.4,"Y");
    TLine *l = new TLine(-2.4,1,2.4,1);
    l->SetLineStyle(2);
    
    l->Draw();   
}
