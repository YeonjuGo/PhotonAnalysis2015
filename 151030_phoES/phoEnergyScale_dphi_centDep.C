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

void phoEnergyScale_dphi_centDep(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer",int ptThr=20, int evpOrder=2){
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

    float epBin[] = {0.0,0.5,0.6,0.7,0.8,0.9,0.95,1.0};
    const int nEP = sizeof(epBin)/sizeof(float)-1;

    int centBin[] = {0,60,200};
    //int centBin[] = {0,20,60,100,200};
    const int nCent = sizeof(centBin)/sizeof(int)-1;
    const char* centSt[nCent];
//    for(int icent=0;icent<nCent;icent++){
//        centSt[icent]=Form("%d-%d%",(int)centBin[icent]/2,(int)centBin[icent+1]/2);
//    }
    centSt[0] = "0-30%";
    centSt[1] = "30-100%";

    TH1D* h1D_ES[nEta][nEP][nCent];//[etaBin][ep][cent]
    TH1D* h1D_ES_ep[nEta][nCent];
    double mean[nEta][nEP][nCent], var[nEta][nEP][nCent], resol[nEta][nEP][nCent], resolVar[nEta][nEP][nCent];

    TCanvas* c2[nEta][nCent];
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c2[ieta][icent] = new TCanvas(Form("c2_ieta%d_icent%d",ieta,icent), Form("pt/refPt distribution %.2f<|eta|<%.2f  %s",etaBin[ieta],etaBin[ieta+1],centSt[icent]), 1200, 600); 
            makeMultiPanelCanvas(c2[ieta][icent],4,2,0.0,0.0,0.2,0.15,0.02);
            h1D_ES_ep[ieta][icent] = new TH1D(Form("h1D_ES_ep_ieta%d_icent%d",ieta,icent),Form(";#Delta#phi(#gamma,#Psi_{%d}); p_{T}^{RECO}/p_{T}^{GEN}",evpOrder), nEP, epBin);
            for(int iep=0; iep<nEP; iep++){
                const char* etaCut = Form("abs(phoEta)>=%.2f && abs(phoEta)<%.2f",etaBin[ieta],etaBin[ieta+1]);
                const char* ptCut = Form("phoEt>=%d",ptThr);
                const char* epCut = Form("abs(phoDPhi_evtpl%d)>=%.2f && abs(phoDPhi_evtpl%d)<%.2f",evpOrder,epBin[iep],evpOrder,epBin[iep+1]);
                const char* centCut = Form("hiBin>=%d && hiBin<%d",centBin[icent],centBin[icent+1]);
                c2[ieta][icent]->cd(iep+1);
                h1D_ES[ieta][iep][icent] = new TH1D(Form("h1D_ES_ieta%d_iep%d_icent%d",ieta,iep,icent),Form("%.1f<#Delta#phi(#gamma,#Psi_{%d})<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",epBin[iep],evpOrder,epBin[iep+1]),40,0,3); 
                tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta][iep][icent]->GetName()), Form("%s && %s && %s && %s",etaCut,ptCut,epCut,centCut));
                h1D_ES[ieta][iep][icent] = (TH1D*)gDirectory->Get(h1D_ES[ieta][iep][icent]->GetName());

                TF1* ff = cleverGaus(h1D_ES[ieta][iep][icent],"h",2.5);
                //gPad->SetLogy();
                mean[ieta][iep][icent] = ff->GetParameter(1);
                var[ieta][iep][icent] = ff->GetParError(1);
                resol[ieta][iep][icent] = ff->GetParameter(2);
                resolVar[ieta][iep][icent] = ff->GetParError(2);
                h1D_ES_ep[ieta][icent]->SetBinContent(iep+1,mean[ieta][iep][icent]);
                h1D_ES_ep[ieta][icent]->SetBinError(iep+1,var[ieta][iep][icent]);

                h1D_ES[ieta][iep][icent]->Draw();
                if(iep==0) {
                    drawText(Form("%s",treePath.Data()),0.6,0.6);
                    drawText(Form("%s",centSt[icent]),0.6,0.68);
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
            if(ieta!=nEta-1) h1D_ES_ep[ieta][icent] -> SetAxisRange(0.9, 1.6, "Y");
            else h1D_ES_ep[ieta][icent] -> SetAxisRange(0.9, 1.6, "Y");
            h1D_ES_ep[ieta][icent] -> SetLineColor(2-icent);
            if(ieta==0) leg->AddEntry(h1D_ES_ep[ieta][icent],centSt[icent]);
            if(icent==0) h1D_ES_ep[ieta][icent] -> Draw();
            else h1D_ES_ep[ieta][icent] -> Draw("same");
            float xpos(0.65),ypos(0.75),dy(0.04);
            if(ieta==0){
                drawText(Form("p_{T}>%d GeV", ptThr),xpos,ypos);
                drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy);
                leg->Draw("same");
            }
            drawText(Form("%.2f<|#eta|<%.2f",etaBin[ieta],etaBin[ieta+1]),xpos,ypos-2*dy); 
            jumSun(epBin[0],1.0,epBin[nEP],1.0);
        }
    }
    c1->SaveAs(Form("png/ES_dphi%d_%s_%s.png",evpOrder,sample.Data(), treePath.Data()));
    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0; icent<nCent; icent++){
            c2[ieta][icent]->SaveAs(Form("png/ES_dphi%dDist_%s_%s_ieta%d_%s.png",evpOrder,sample.Data(),treePath.Data(),ieta,centSt[icent]));
        }
    }
}
