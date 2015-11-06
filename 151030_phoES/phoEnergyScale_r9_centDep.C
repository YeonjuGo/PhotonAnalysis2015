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

void phoEnergyScale_r9_centDep(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer",int ptThr=20){
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

    float r9Bin[] = {0.0,0.5,0.6,0.7,0.8,0.9,0.95,1.0};
    const int nR9 = sizeof(r9Bin)/sizeof(float)-1;

    int centBin[] = {0,60,200};
    //int centBin[] = {0,20,60,100,200};
    const int nCent = sizeof(centBin)/sizeof(int)-1;
    const char* centSt[nCent];
//    for(int icent=0;icent<nCent;icent++){
//        centSt[icent]=Form("%d-%d%",(int)centBin[icent]/2,(int)centBin[icent+1]/2);
//    }
    centSt[0] = "0-30%";
    centSt[1] = "30-100%";

    TH1D* h1D_ES[nEta][nR9][nCent];//[etaBin][R9][cent]
    TH1D* h1D_ES_r9[nEta][nCent];
    double mean[nEta][nR9][nCent], var[nEta][nR9][nCent], resol[nEta][nR9][nCent], resolVar[nEta][nR9][nCent];

    TCanvas* c2[nEta][nCent];
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c2[ieta][icent] = new TCanvas(Form("c2_ieta%d_icent%d",ieta,icent), Form("pt/refPt distribution %.2f<|eta|<%.2f  %s",etaBin[ieta],etaBin[ieta+1],centSt[icent]), 1200, 600); 
            makeMultiPanelCanvas(c2[ieta][icent],4,2,0.0,0.0,0.2,0.15,0.02);
            h1D_ES_r9[ieta][icent] = new TH1D(Form("h1D_ES_r9_ieta%d_icent%d",ieta,icent),";photon R9; p_{T}^{RECO}/p_{T}^{GEN}", nR9, r9Bin);
            for(int ir9=0; ir9<nR9; ir9++){
                const char* etaCut = Form("abs(phoEta)>=%.2f && abs(phoEta)<%.2f",etaBin[ieta],etaBin[ieta+1]);
                const char* ptCut = Form("phoEt>=%d",ptThr);
                const char* r9Cut = Form("phoR9>=%.2f && phoR9<%.2f",r9Bin[ir9],r9Bin[ir9+1]);
                const char* centCut = Form("hiBin>=%d && hiBin<%d",centBin[icent],centBin[icent+1]);
                c2[ieta][icent]->cd(ir9+1);
                h1D_ES[ieta][ir9][icent] = new TH1D(Form("h1D_ES_ieta%d_ir9%d_icent%d",ieta,ir9,icent),Form("%.1f<photon R9<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",r9Bin[ir9],r9Bin[ir9+1]),40,0,2); 
                tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta][ir9][icent]->GetName()), Form("%s && %s && %s && %s",etaCut,ptCut,r9Cut,centCut));
                h1D_ES[ieta][ir9][icent] = (TH1D*)gDirectory->Get(h1D_ES[ieta][ir9][icent]->GetName());

                TF1* ff = cleverGaus(h1D_ES[ieta][ir9][icent],"h",2.5);
                //gPad->SetLogy();
                mean[ieta][ir9][icent] = ff->GetParameter(1);
                var[ieta][ir9][icent] = ff->GetParError(1);
                resol[ieta][ir9][icent] = ff->GetParameter(2);
                resolVar[ieta][ir9][icent] = ff->GetParError(2);
                h1D_ES_r9[ieta][icent]->SetBinContent(ir9+1,mean[ieta][ir9][icent]);
                h1D_ES_r9[ieta][icent]->SetBinError(ir9+1,var[ieta][ir9][icent]);

                h1D_ES[ieta][ir9][icent]->Draw();
                if(ir9==0) {
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
            if(ieta!=nEta-1) h1D_ES_r9[ieta][icent] -> SetAxisRange(0.9, 1.6, "Y");
            else h1D_ES_r9[ieta][icent] -> SetAxisRange(0.9, 1.6, "Y");
            h1D_ES_r9[ieta][icent] -> SetLineColor(2-icent);
            if(ieta==0) leg->AddEntry(h1D_ES_r9[ieta][icent],centSt[icent]);
            if(icent==0) h1D_ES_r9[ieta][icent] -> Draw();
            else h1D_ES_r9[ieta][icent] -> Draw("same");
            float xpos(0.65),ypos(0.75),dy(0.04);
            if(ieta==0){
                drawText(Form("p_{T}>%d GeV", ptThr),xpos,ypos);
                drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy);
                leg->Draw("same");
            }
            drawText(Form("%.2f<|#eta|<%.2f",etaBin[ieta],etaBin[ieta+1]),xpos,ypos-2*dy); 
            jumSun(r9Bin[0],1.0,r9Bin[nR9],1.0);
        }
    }
    c1->SaveAs(Form("png/ES_r9_%s_%s.png",sample.Data(), treePath.Data()));
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c2[ieta][icent]->SaveAs(Form("png/ES_r9Dist_%s_%s_ieta%d_%s.png",sample.Data(),treePath.Data(),ieta,centSt[icent]));
        }
    }
}
