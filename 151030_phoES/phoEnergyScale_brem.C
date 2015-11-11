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

void phoEnergyScale_brem(TString sample = "SingleGammaFlatPt10To200", TString treePath="ggHiNtuplizerGED",int ptThr=20){

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
    TString inputFileName = Form("/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/150908_genMatchingSkim/skimFiles/jskim_%s_%s_genMatched_inclusivePho.root",sample.Data(),treePath.Data());
    TFile* fin = new TFile(inputFileName);
    TTree* tr = (TTree*)fin->Get("t_pho");

    float etaBin[] = {0.00,1.44,2.00,3.00};
    //float etaBin[] = {0.00,1.44,2.00,2.40,3.00};
    const int nEta = sizeof(etaBin)/sizeof(float)-1;

    //brem = phiWidth/etaWidth
    float bremBin[] = {0.0,1.0,1.5,2.0,3.0,16.0};
    const int nBrem = sizeof(bremBin)/sizeof(float)-1;
    TH1D* h1D_ES[nEta][nBrem];//[etaBin][eventPlane]
    TH1D* h1D_ES_brem[nEta];
    double mean[nEta][nBrem], var[nEta][nBrem], resol[nEta][nBrem], resolVar[nEta][nBrem];
    
    TCanvas* c2[nEta];
    for(int ieta=0; ieta<nEta; ieta++){
        c2[ieta] = new TCanvas(Form("c2_ieta%d",ieta), Form("pt/refPt distribution %.2f<|eta|<%.2f",etaBin[ieta],etaBin[ieta+1]), 1200, 600); 
        makeMultiPanelCanvas(c2[ieta],4,2,0.0,0.0,0.2,0.15,0.02);
        h1D_ES_brem[ieta] = new TH1D(Form("h1D_ES_brem_ieta%d",ieta),";photon brem; p_{T}^{RECO}/p_{T}^{GEN}", nBrem, bremBin);
        for(int ibrem=0; ibrem<nBrem; ibrem++){
            const char* etaCut = Form("abs(phoEta)>=%.2f && abs(phoEta)<%.2f",etaBin[ieta],etaBin[ieta+1]);
            const char* ptCut = Form("phoEt>=%d",ptThr);
            const char* bremCut = Form("phoSCBrem>=%.2f && phoSCBrem<%.2f",bremBin[ibrem],bremBin[ibrem+1]);
            c2[ieta]->cd(ibrem+1);
            h1D_ES[ieta][ibrem] = new TH1D(Form("h1D_ES_ieta%d_ibrem%d",ieta,ibrem),Form("%.1f<photon brem<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",bremBin[ibrem],bremBin[ibrem+1]),40,0,3); 
            tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta][ibrem]->GetName()), Form("%s && %s && %s",etaCut,ptCut,bremCut));
            h1D_ES[ieta][ibrem] = (TH1D*)gDirectory->Get(h1D_ES[ieta][ibrem]->GetName());

            TF1* ff = cleverGaus(h1D_ES[ieta][ibrem],"h",2.5);
            //gPad->SetLogy();
            mean[ieta][ibrem] = ff->GetParameter(1);
            var[ieta][ibrem] = ff->GetParError(1);
            resol[ieta][ibrem] = ff->GetParameter(2);
            resolVar[ieta][ibrem] = ff->GetParError(2);
            h1D_ES_brem[ieta]->SetBinContent(ibrem+1,mean[ieta][ibrem]);
            h1D_ES_brem[ieta]->SetBinError(ibrem+1,var[ieta][ibrem]);

            h1D_ES[ieta][ibrem]->Draw();
            if(ibrem==0) drawText(Form("%s",treePath.Data()),0.6,0.6);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Energy Scale vs. pt", 1000,600); 
    c1->Divide(2,2);
    for(int ieta=0; ieta<nEta; ieta++){
        c1->cd(ieta+1);
        if(ieta!=nEta-1) h1D_ES_brem[ieta] -> SetAxisRange(0.9, 1.3, "Y");
        else h1D_ES_brem[ieta] -> SetAxisRange(0.8, 1.3, "Y");
        h1D_ES_brem[ieta] -> Draw();
        float xpos(0.65),ypos(0.75),dy(0.04);
        if(ieta==0){
            drawText(Form("p_{T}>%d GeV", ptThr),xpos,ypos);
            drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy);
        }
        drawText(Form("%.2f<|#eta|<%.2f",etaBin[ieta],etaBin[ieta+1]),xpos,ypos-2*dy); 
        jumSun(bremBin[0],1.0,bremBin[nBrem],1.0);
    }
    c1->SaveAs(Form("png/ES_brem_%s_%s.png",sample.Data(), treePath.Data()));
    for(int ieta=0; ieta<nEta; ieta++){
        c2[ieta]->SaveAs(Form("png/ES_bremDist_%s_%s_ieta%d.png",sample.Data(),treePath.Data(),ieta));
    }
}
