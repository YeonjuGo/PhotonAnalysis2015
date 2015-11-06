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

void phoEnergyScale_pt(TString sample = "SingleGammaFlatPt10To200", TString treePath="ggHiNtuplizerGED"){
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

    float ptBin[] = {20.0,30.0,40.0,50.0,70.0,90.0,110.0};
    const int nPt = sizeof(ptBin)/sizeof(float)-1;
    TH1D* h1D_ES[nEta][nPt];//[etaBin][pt]
    TH1D* h1D_ES_pt[nEta];
    double mean[nEta][nPt], var[nEta][nPt], resol[nEta][nPt], resolVar[nEta][nPt];
    
    TCanvas* c2[nEta];
    for(int ieta=0; ieta<nEta; ieta++){
        c2[ieta] = new TCanvas(Form("c2_ieta%d",ieta), Form("pt/refPt distribution %.2f<|eta|<%.2f",etaBin[ieta],etaBin[ieta+1]), 1200, 600); 
        makeMultiPanelCanvas(c2[ieta],4,2,0.0,0.0,0.2,0.15,0.02);
        h1D_ES_pt[ieta] = new TH1D(Form("h1D_ES_pt_ieta%d",ieta),";photon p_{T}^{GEN}; p_{T}^{RECO}/p_{T}^{GEN}", nPt, ptBin);
        for(int ipt=0; ipt<nPt; ipt++){
            c2[ieta]->cd(ipt+1);
            h1D_ES[ieta][ipt] = new TH1D(Form("h1D_ES_ieta%d_ipt%d",ieta,ipt),Form("%.1f<p_{T}^{GEN}<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",ptBin[ipt],ptBin[ipt+1]),40,0,3); 
            tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta][ipt]->GetName()), Form("abs(phoEta)>=%.2f && abs(phoEta)<%.2f && phoEt>=%.1f && phoEt<%.1f",etaBin[ieta],etaBin[ieta+1],ptBin[ipt],ptBin[ipt+1]));
            h1D_ES[ieta][ipt] = (TH1D*)gDirectory->Get(h1D_ES[ieta][ipt]->GetName());

            TF1* ff = cleverGaus(h1D_ES[ieta][ipt],"h",2.5);
            //gPad->SetLogy();
            mean[ieta][ipt] = ff->GetParameter(1);
            var[ieta][ipt] = ff->GetParError(1);
            resol[ieta][ipt] = ff->GetParameter(2);
            resolVar[ieta][ipt] = ff->GetParError(2);
            h1D_ES_pt[ieta]->SetBinContent(ipt+1,mean[ieta][ipt]);
            h1D_ES_pt[ieta]->SetBinError(ipt+1,var[ieta][ipt]);

            h1D_ES[ieta][ipt]->Draw();
            if(ipt==0) drawText(Form("%s",treePath.Data()),0.6,0.6);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Energy Scale vs. pt", 1000,600); 
    c1->Divide(2,2);
    for(int ieta=0; ieta<nEta; ieta++){
        c1->cd(ieta+1);
        if(ieta!=4) h1D_ES_pt[ieta] -> SetAxisRange(0.9, 1.3, "Y");
        else h1D_ES_pt[ieta] -> SetAxisRange(0.9, 3.0, "Y");
        h1D_ES_pt[ieta] -> Draw();
        float xpos(0.65),ypos(0.75),dy(0.04);
        if(ieta==0){
            //drawText(Form("p_{T}>%d GeV", ptThr),xpos,ypos);
            drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy);
        }
        drawText(Form("%.2f<|#eta|<%.2f",etaBin[ieta],etaBin[ieta+1]),xpos,ypos-2*dy); 
        jumSun(20.0,1.0,ptBin[nPt],1.0);
    }
    c1->SaveAs(Form("png/ES_pt_%s_%s.png",sample.Data(), treePath.Data()));
    for(int ieta=0; ieta<nEta; ieta++){
        c2[ieta]->SaveAs(Form("png/ES_ptDist_%s_%s_ieta%d.png",sample.Data(),treePath.Data(),ieta));
    }
}
