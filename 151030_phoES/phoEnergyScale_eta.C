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

void phoEnergyScale_eta(TString sample = "SingleGammaFlatPt10To200", TString treePath="ggHiNtuplizerGED"){
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

    double etaBin[] = {-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0};
    const int nEta = sizeof(etaBin)/sizeof(double)-1;

    double ptThrBin[] = {20.0,40.0,60.0};
    const int nPtThr = sizeof(ptThrBin)/sizeof(double);
    TH1D* h1D_ES[nPtThr][nEta];//[ptThr][etaBin]
    TH1D* h1D_ES_eta[nPtThr];
    TCanvas* c2[nPtThr];
    double mean[nPtThr][nEta], var[nPtThr][nEta], resol[nPtThr][nEta], resolVar[nPtThr][nEta];

    for(int iptThr=0; iptThr<nPtThr; iptThr++){
        h1D_ES_eta[iptThr] = new TH1D(Form("h1D_ES_iptThr%d",iptThr),";photon #eta; p_{T}^{RECO}/p_{T}^{GEN}", nEta, etaBin);
        c2[iptThr] = new TCanvas(Form("c2_iptThr%d",iptThr), Form("recoPt/genPt distribution pT>%d GeV",(int)ptThrBin[iptThr]), 1200, 900); 
        makeMultiPanelCanvas(c2[iptThr],4,3,0.0,0.0,0.2,0.15,0.02);

        for(int ieta=0; ieta<nEta; ieta++){
            c2[iptThr]->cd(ieta+1);
            h1D_ES[iptThr][ieta] = new TH1D(Form("h1D_ES_ieta%d_iptThr%d",ieta,iptThr),Form("%.1f<|#eta|<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",etaBin[ieta],etaBin[ieta+1]),40,0,3); 
            tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[iptThr][ieta]->GetName()), Form("phoEta>=%.1f && phoEta<%.1f && phoEt>=%d",etaBin[ieta],etaBin[ieta+1],(int)ptThrBin[iptThr]));
            h1D_ES[iptThr][ieta] = (TH1D*)gDirectory->Get(h1D_ES[iptThr][ieta]->GetName());

            TF1* ff = cleverGaus(h1D_ES[iptThr][ieta],"h",2.5);
            //gPad->SetLogy();
            mean[iptThr][ieta] = ff->GetParameter(1);
            var[iptThr][ieta] = ff->GetParError(1);
            resol[iptThr][ieta] = ff->GetParameter(2);
            resolVar[iptThr][ieta] = ff->GetParError(2);
            h1D_ES_eta[iptThr]->SetBinContent(ieta+1,mean[iptThr][ieta]);
            h1D_ES_eta[iptThr]->SetBinError(ieta+1,var[iptThr][ieta]);

            h1D_ES[iptThr][ieta]->Draw();
            if(ieta==0) drawText(Form("%s",treePath.Data()),0.2,0.85);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Energy Scale vs. photon eta", 1000,600); 
    c1->Divide(2,2);
    for(int iptThr=0; iptThr<nPtThr; iptThr++){
        c1->cd(iptThr+1);
        h1D_ES_eta[iptThr] -> SetAxisRange(0.9, 1.3, "Y");
        h1D_ES_eta[iptThr] -> Draw();
        float xpos(0.65),ypos(0.75),dy(0.04);
        if(iptThr==0){
            drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy); 
        }
        drawText(Form("p_{T}>%d GeV", (int)ptThrBin[iptThr]),xpos,ypos);
        jumSun(etaBin[0],1,etaBin[nEta],1);
    }
    c1->SaveAs(Form("png/ES_eta_%s_%s.png",sample.Data(),treePath.Data()));
    for(int iptThr=0; iptThr<nPtThr; iptThr++){
        c2[iptThr]->SaveAs(Form("png/ES_etaDist_%s_%s_ptThr%d.png",sample.Data(),treePath.Data(),(int)ptThrBin[iptThr]));
    }
}
