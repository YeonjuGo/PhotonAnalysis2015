// yeonju
// pt fine bin 
//#include "../HiForestAnalysis/hiForest.h"
//#include "../gammaJetAnalysis/CutAndBinCollection2012.h"
#include <TStyle.h>
#include <TH1D.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <iostream>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TCut.h>
#include "../yjUtility.h"

void phoEnergyScale_eta(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer", int ptThr=20){
    TH1::SetDefaultSumw2();
    // gStyle->SetOptFit(0);
    gStyle -> SetOptStat(0);
    //  gStyle -> SetTitleYOffset(2.35);
    gStyle -> SetTitleYSize(0.04);
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");

   TString inputFileName = Form("/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/150908_voronoiStudy/skimFiles/jskim_%s_%s_genMatched_inclusivePho.root",sample.Data(),treePath.Data());
    TFile* fin = new TFile(inputFileName);
    TTree* tr = (TTree*)fin->Get("t_pho");

    double etaBin[] = {-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0};
    const int nEta = sizeof(etaBin)/sizeof(double)-1;

    TH1D* h1D_ES[nEta];//[etaBin]
    TH1D* h1D_ES_eta = new TH1D(Form("h1D_ES_ptThr%d",ptThr),";photon #eta; p_{T}^{RECO}/p_{T}^{GEN}", nEta, etaBin);
    TCanvas* c2 = new TCanvas("c2", "pt/refPt distribution", 1200, 900); 
    makeMultiPanelCanvas(c2,4,3,0.0,0.0,0.2,0.15,0.02);
    double mean[nEta], var[nEta], resol[nEta], resolVar[nEta];

    for(int ieta=0; ieta<nEta; ieta++){
        c2->cd(ieta+1);
        h1D_ES[ieta] = new TH1D(Form("h1D_ES_ieta%d",ieta),Form("%.1f<|#eta|<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",etaBin[ieta],etaBin[ieta+1]),40,0,3); 
        tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta]->GetName()), Form("phoEta>=%.1f && phoEta<%.1f && phoEt>=%d",etaBin[ieta],etaBin[ieta+1],ptThr));
        h1D_ES[ieta] = (TH1D*)gDirectory->Get(h1D_ES[ieta]->GetName());

        TF1* ff = cleverGaus(h1D_ES[ieta],"h",1.5);
        //gPad->SetLogy();
        mean[ieta] = ff->GetParameter(1);
        var[ieta] = ff->GetParError(1);
        resol[ieta] = ff->GetParameter(2);
        resolVar[ieta] = ff->GetParError(2);
        h1D_ES_eta->SetBinContent(ieta+1,mean[ieta]);
        h1D_ES_eta->SetBinError(ieta+1,var[ieta]);

        h1D_ES[ieta]->Draw();
        if(ieta==0) {
            drawText(Form("%s",treePath.Data()),0.2,0.85);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Energy Scale vs. photon eta", 600,300); 
    h1D_ES_eta -> SetAxisRange(0.9, 1.3, "Y");
    h1D_ES_eta -> Draw();
    float xpos(0.6),ypos(0.8),dy(0.06);
    drawText(Form("p_{T}>%d GeV", ptThr),xpos,ypos);
    drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy); 
    jumSun(-3.0,1,3.0,1);

    c1->SaveAs(Form("png/ES_eta_%s_%s_ptThr%d.png",sample.Data(),treePath.Data(),ptThr));
    c2->SaveAs(Form("png/ES_etaDist_%s_%s_ptThr%d.png",sample.Data(),treePath.Data(),ptThr));
}
