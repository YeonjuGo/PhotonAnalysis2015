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

void phoEnergyScale_pt_centDep(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizerGED"){
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

    int centBin[] = {0,60,200};
    const int nCent = sizeof(centBin)/sizeof(int)-1;
    const char* centSt[nCent];
    centSt[0] = "0-30%";
    centSt[1] = "30-100%";

    TH1D* h1D_ES[nEta][nPt][nCent];//[etaBin][pt][cent]
    TH1D* h1D_ES_pt[nEta][nCent];
    double mean[nEta][nPt][nCent], var[nEta][nPt][nCent], resol[nEta][nPt][nCent], resolVar[nEta][nPt][nCent];

    TCanvas* c2[nEta][nCent];
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c2[ieta][icent] = new TCanvas(Form("c2_ieta%d_icent%d",ieta,icent), Form("pt/refPt distribution %.2f<|eta|<%.2f  %s",etaBin[ieta],etaBin[ieta+1],centSt[icent]), 1200, 600); 
            makeMultiPanelCanvas(c2[ieta][icent],4,2,0.0,0.0,0.2,0.15,0.02);
            h1D_ES_pt[ieta][icent] = new TH1D(Form("h1D_ES_pt_ieta%d_icent%d",ieta,icent),";photon pt; p_{T}^{RECO}/p_{T}^{GEN}", nPt, ptBin);
            for(int ipt=0; ipt<nPt; ipt++){
                const char* etaCut = Form("abs(phoEta)>=%.2f && abs(phoEta)<%.2f",etaBin[ieta],etaBin[ieta+1]);
                const char* ptCut = Form("phoEt>=%.2f && phoEt<%.2f",ptBin[ipt],ptBin[ipt+1]);
                const char* centCut = Form("hiBin>=%d && hiBin<%d",centBin[icent],centBin[icent+1]);
                c2[ieta][icent]->cd(ipt+1);
                h1D_ES[ieta][ipt][icent] = new TH1D(Form("h1D_ES_ieta%d_ipt%d_icent%d",ieta,ipt,icent),Form("%.1f<photon pt<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",ptBin[ipt],ptBin[ipt+1]),40,0,3); 
                tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta][ipt][icent]->GetName()), Form("%s && %s && %s && %s",etaCut,ptCut,ptCut,centCut));
                h1D_ES[ieta][ipt][icent] = (TH1D*)gDirectory->Get(h1D_ES[ieta][ipt][icent]->GetName());

                TF1* ff = cleverGaus(h1D_ES[ieta][ipt][icent],"h",2.5);
                //gPad->SetLogy();
                mean[ieta][ipt][icent] = ff->GetParameter(1);
                var[ieta][ipt][icent] = ff->GetParError(1);
                resol[ieta][ipt][icent] = ff->GetParameter(2);
                resolVar[ieta][ipt][icent] = ff->GetParError(2);
                h1D_ES_pt[ieta][icent]->SetBinContent(ipt+1,mean[ieta][ipt][icent]);
                h1D_ES_pt[ieta][icent]->SetBinError(ipt+1,var[ieta][ipt][icent]);

                h1D_ES[ieta][ipt][icent]->Draw();
                if(ipt==0) {
                    drawText(Form("%s",treePath.Data()),0.6,0.6);
                    drawText(Form("%s",centSt[icent]),0.2,0.78);
                }
            }
        }
    }

    TLegend* leg = new TLegend(0.3929766,0.6934307,0.6923077,0.9927007,NULL,"brNDC");
    legStyle(leg);

    TCanvas* c1 = new TCanvas("c1", "Energy Scale vs. pt", 1000,600); 
    c1->Divide(2,2);
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c1->cd(ieta+1);
            if(ieta!=nEta-1) h1D_ES_pt[ieta][icent] -> SetAxisRange(0.9, 1.5, "Y");
            else h1D_ES_pt[ieta][icent] -> SetAxisRange(0.9, 1.5, "Y");
            h1D_ES_pt[ieta][icent] -> SetLineColor(2-icent);
            if(ieta==0) leg->AddEntry(h1D_ES_pt[ieta][icent],centSt[icent]);
            if(icent==0) h1D_ES_pt[ieta][icent] -> Draw();
            else h1D_ES_pt[ieta][icent] -> Draw("same");
            float xpos(0.65),ypos(0.75),dy(0.04);
            if(ieta==0){
                drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy);
                leg->Draw("same");
            }
            drawText(Form("%.2f<|#eta|<%.2f",etaBin[ieta],etaBin[ieta+1]),xpos,ypos-2*dy); 
            jumSun(ptBin[0],1.0,ptBin[nPt],1.0);
        }
    }
    c1->SaveAs(Form("png/ES_pt_%s_%s.png",sample.Data(), treePath.Data()));
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c2[ieta][icent]->SaveAs(Form("png/ES_ptDist_%s_%s_ieta%d_%s.png",sample.Data(),treePath.Data(),ieta,centSt[icent]));
        }
    }
}
