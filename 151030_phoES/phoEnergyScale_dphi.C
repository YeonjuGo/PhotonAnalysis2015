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

void phoEnergyScale_dphi(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer", int ptThr=20, int evpOrder=2){
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

    double etaBin[] = {0.00,1.44,2.00,3.00};
    //double etaBin[] = {0.00,1.44,2.00,2.40,3.00};
    const int nEta = sizeof(etaBin)/sizeof(double)-1;

    const int nEP = 8;
    const double epmax = TMath::Pi();
    double epBin[nEP+1];
    for(int i=0; i<nEP+1;i++){
        double binWidth = (epmax)/nEP;
        epBin[i] = 0.0 + binWidth*i;
    }

    TH1D* h1D_ES[nEta][nEP];//[etaBin][eventPlane]
    TH1D* h1D_ES_dphi[nEta];
    double mean[nEta][nEP], var[nEta][nEP], resol[nEta][nEP], resolVar[nEta][nEP];
    
    TCanvas* c2[nEta];
    for(int ieta=0; ieta<nEta; ieta++){
        c2[ieta] = new TCanvas(Form("c2_ieta%d",ieta), "pt/refPt distribution", 1200, 600); 
        makeMultiPanelCanvas(c2[ieta],4,2,0.0,0.0,0.2,0.15,0.02);
        h1D_ES_dphi[ieta] = new TH1D(Form("h1D_ES_dphi_ptThr%d_evpOrder%d_ieta%d",ptThr,evpOrder,ieta),";#Delta#phi(#gamma,#Psi_{2}); p_{T}^{RECO}/p_{T}^{GEN}", nEP, epBin);
        for(int iep=0; iep<nEP; iep++){
            c2[ieta]->cd(iep+1);
            h1D_ES[ieta][iep] = new TH1D(Form("h1D_ES_ieta%d_iep%d",ieta,iep),Form("%.1f<#Delta#phi<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",epBin[iep],epBin[iep+1]),40,0,3); 
            tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta][iep]->GetName()), Form("phoEta>=%.2f && phoEta<%.2f && phoEt>=%d && abs(phoDPhi_evtpl%d)>=%.2f && abs(phoDPhi_evtpl%d)<%.2f",etaBin[ieta],etaBin[ieta+1],ptThr,evpOrder,epBin[iep],evpOrder,epBin[iep+1]));
            h1D_ES[ieta][iep] = (TH1D*)gDirectory->Get(h1D_ES[ieta][iep]->GetName());

            TF1* ff = cleverGaus(h1D_ES[ieta][iep],"h",1.5);
            //gPad->SetLogy();
            mean[ieta][iep] = ff->GetParameter(1);
            var[ieta][iep] = ff->GetParError(1);
            resol[ieta][iep] = ff->GetParameter(2);
            resolVar[ieta][iep] = ff->GetParError(2);
            h1D_ES_dphi[ieta]->SetBinContent(iep+1,mean[ieta][iep]);
            h1D_ES_dphi[ieta]->SetBinError(iep+1,var[ieta][iep]);

            h1D_ES[ieta][iep]->Draw();
            if(ieta==0) drawText(Form("%s",treePath.Data()),0.6,0.6);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Energy Scale vs. dphi", 1000,600); 
    c1->Divide(2,2);
    for(int ieta=0; ieta<nEta; ieta++){
        c1->cd(ieta+1);
        if(ieta!=4) h1D_ES_dphi[ieta] -> SetAxisRange(0.9, 1.3, "Y");
        else h1D_ES_dphi[ieta] -> SetAxisRange(0.9, 3.0, "Y");
        h1D_ES_dphi[ieta] -> Draw();
        float xpos(0.6),ypos(0.6),dy(0.06);
        if(ieta==0){
            drawText(Form("p_{T}>%d GeV", ptThr),xpos,ypos);
            drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy);
        }
        drawText(Form("%.2f<|#eta|<%.2f",etaBin[ieta],etaBin[ieta+1]),xpos,ypos-2*dy); 
        jumSun(0.0,1.0,epmax,1.0);
    }
    c1->SaveAs(Form("png/ES_dphi%d_%s_%s_ptThr%d.png",evpOrder,sample.Data(),treePath.Data(),ptThr));
    for(int ieta=0; ieta<nEta; ieta++){
        c2[ieta]->SaveAs(Form("png/ES_dphi%dDist_%s_%s_ptThr%d_ieta%d.png",evpOrder,sample.Data(),treePath.Data(),ptThr,ieta));
    }
}
