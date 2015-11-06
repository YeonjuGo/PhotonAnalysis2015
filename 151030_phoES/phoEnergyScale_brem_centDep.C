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

void phoEnergyScale_brem_centDep(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizerGED",int ptThr=20){
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

    float bremBin[] = {0.0,1.0,1.5,2.0,3.0,16.0}; 
    const int nBrem = sizeof(bremBin)/sizeof(float)-1;

    int centBin[] = {0,60,200};
    const int nCent = sizeof(centBin)/sizeof(int)-1;
    const char* centSt[nCent];
    centSt[0] = "0-30%";
    centSt[1] = "30-100%";

    TH1D* h1D_ES[nEta][nBrem][nCent];//[etaBin][brem][cent]
    TH1D* h1D_ES_brem[nEta][nCent];
    double mean[nEta][nBrem][nCent], var[nEta][nBrem][nCent], resol[nEta][nBrem][nCent], resolVar[nEta][nBrem][nCent];

    TCanvas* c2[nEta][nCent];
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c2[ieta][icent] = new TCanvas(Form("c2_ieta%d_icent%d",ieta,icent), Form("pt/refPt distribution %.2f<|eta|<%.2f  %s",etaBin[ieta],etaBin[ieta+1],centSt[icent]), 1200, 600); 
            makeMultiPanelCanvas(c2[ieta][icent],4,2,0.0,0.0,0.2,0.15,0.02);
            h1D_ES_brem[ieta][icent] = new TH1D(Form("h1D_ES_brem_ieta%d_icent%d",ieta,icent),";photon brem; p_{T}^{RECO}/p_{T}^{GEN}", nBrem, bremBin);
            for(int ibrem=0; ibrem<nBrem; ibrem++){
                const char* etaCut = Form("abs(phoEta)>=%.2f && abs(phoEta)<%.2f",etaBin[ieta],etaBin[ieta+1]);
                const char* ptCut = Form("phoEt>=%d",ptThr);
                const char* bremCut = Form("phoSCBrem>=%.2f && phoSCBrem<%.2f",bremBin[ibrem],bremBin[ibrem+1]);
                const char* centCut = Form("hiBin>=%d && hiBin<%d",centBin[icent],centBin[icent+1]);
                c2[ieta][icent]->cd(ibrem+1);
                h1D_ES[ieta][ibrem][icent] = new TH1D(Form("h1D_ES_ieta%d_ibrem%d_icent%d",ieta,ibrem,icent),Form("%.1f<photon SC Brem(=PhiWidth/EtaWidth)<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",bremBin[ibrem],bremBin[ibrem+1]),40,0,3); 
                tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta][ibrem][icent]->GetName()), Form("%s && %s && %s && %s",etaCut,ptCut,bremCut,centCut));
                h1D_ES[ieta][ibrem][icent] = (TH1D*)gDirectory->Get(h1D_ES[ieta][ibrem][icent]->GetName());

                TF1* ff = cleverGaus(h1D_ES[ieta][ibrem][icent],"h",1.5);
                //gPad->SetLogy();
                mean[ieta][ibrem][icent] = ff->GetParameter(1);
                var[ieta][ibrem][icent] = ff->GetParError(1);
                resol[ieta][ibrem][icent] = ff->GetParameter(2);
                resolVar[ieta][ibrem][icent] = ff->GetParError(2);
                h1D_ES_brem[ieta][icent]->SetBinContent(ibrem+1,mean[ieta][ibrem][icent]);
                h1D_ES_brem[ieta][icent]->SetBinError(ibrem+1,var[ieta][ibrem][icent]);

                h1D_ES[ieta][ibrem][icent]->Draw();
                if(ibrem==0) {
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
            if(ieta!=nEta-1) h1D_ES_brem[ieta][icent] -> SetAxisRange(0.9, 1.3, "Y");
            else h1D_ES_brem[ieta][icent] -> SetAxisRange(0.8, 2.0, "Y");
            h1D_ES_brem[ieta][icent] -> SetLineColor(2-icent);
            if(ieta==0) leg->AddEntry(h1D_ES_brem[ieta][icent],centSt[icent]);
            if(icent==0) h1D_ES_brem[ieta][icent] -> Draw();
            else h1D_ES_brem[ieta][icent] -> Draw("same");
            float xpos(0.65),ypos(0.75),dy(0.04);
            if(ieta==0){
                drawText(Form("p_{T}>%d GeV", ptThr),xpos,ypos);
                drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy);
                leg->Draw("same");
            }
            drawText(Form("%.2f<|#eta|<%.2f",etaBin[ieta],etaBin[ieta+1]),xpos,ypos-2*dy); 
            jumSun(bremBin[0],1.0,bremBin[nBrem],1.0);
        }
    }
    c1->SaveAs(Form("png/ES_brem_%s_%s.png",sample.Data(), treePath.Data()));
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c2[ieta][icent]->SaveAs(Form("png/ES_bremDist_%s_%s_ieta%d_%s.png",sample.Data(),treePath.Data(),ieta,centSt[icent]));
        }
    }
}
