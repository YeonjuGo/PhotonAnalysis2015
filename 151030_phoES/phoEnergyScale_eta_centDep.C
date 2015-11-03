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

void phoEnergyScale_eta_centDep(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer", int ptThr=20){
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

    int centBin[] = {0,60,100};
    const int nCent = sizeof(centBin)/sizeof(int)-1;
    const char* centSt[nCent];
    centSt[0] = "0-30%";
    centSt[1] = "30-100%";    

    TH1D* h1D_ES[nEta][nCent];//[etaBin][centBin]
    TH1D* h1D_ES_eta[nCent];
    TCanvas* c2[nCent];
    for(int icent=0;icent<nCent;icent++){ 
        h1D_ES_eta[icent]= new TH1D(Form("h1D_ES_ptThr%d_icent%d",ptThr,icent),";photon #eta; p_{T}^{RECO}/p_{T}^{GEN}", nEta, etaBin);
        c2[icent] = new TCanvas(Form("c2_icent%d",icent), "pt/refPt distribution", 1200, 900); 
        makeMultiPanelCanvas(c2[icent],4,3,0.0,0.0,0.2,0.15,0.02);
    }
    double mean[nEta][nCent], var[nEta][nCent], resol[nEta][nCent], resolVar[nEta][nCent];


    for(int ieta=0; ieta<nEta; ieta++){
        for(int icent=0;icent<nCent;icent++){ 
            c2[icent]->cd(ieta+1);
            const char* centCut = Form("hiBin>=%d && hiBin<%d",centBin[icent],centBin[icent+1]);
            h1D_ES[ieta][icent] = new TH1D(Form("h1D_ES_ieta%d_icent%d",ieta,icent),Form("%.1f<|#eta|<%.1f;p_{T}^{RECO}/p_{T}^{GEN};entries",etaBin[ieta],etaBin[ieta+1]),40,0,3); 
            tr->Draw(Form("phoEt/mcPt>>+%s",h1D_ES[ieta][icent]->GetName()), Form("phoEta>=%.1f && phoEta<%.1f && phoEt>=%d && %s",etaBin[ieta],etaBin[ieta+1],ptThr,centCut));
            h1D_ES[ieta][icent] = (TH1D*)gDirectory->Get(h1D_ES[ieta][icent]->GetName());

            TF1* ff = cleverGaus(h1D_ES[ieta][icent],"h",1.0);
            //gPad->SetLogy();
            mean[ieta][icent] = ff->GetParameter(1);
            var[ieta][icent] = ff->GetParError(1);
            resol[ieta][icent] = ff->GetParameter(2);
            resolVar[ieta][icent] = ff->GetParError(2);
            h1D_ES_eta[icent]->SetBinContent(ieta+1,mean[ieta][icent]);
            h1D_ES_eta[icent]->SetBinError(ieta+1,var[ieta][icent]);

            h1D_ES[ieta][icent]->Draw();
            if(ieta==0) {
                drawText(Form("%s",treePath.Data()),0.2,0.85);
                drawText(Form("%s",centSt[icent]),0.2,0.78);
            }
        }
    }

    TLegend* leg = new TLegend(0.3929766,0.6934307,0.6923077,0.9927007,NULL,"brNDC");
    legStyle(leg);
    TCanvas* c1 = new TCanvas("c1", "Energy Scale vs. photon eta", 600,300); 
    for(int icent=0;icent<nCent;icent++){ 
        h1D_ES_eta[icent] -> SetAxisRange(0.9, 1.3, "Y");
        h1D_ES_eta[icent] -> SetLineColor(2-icent);
        leg->AddEntry(h1D_ES_eta[icent],centSt[icent]);
        if(icent==0) h1D_ES_eta[icent] -> Draw();
        else h1D_ES_eta[icent] -> Draw("same");
    }
    leg->Draw("same");
    float xpos(0.6),ypos(0.8),dy(0.06);
    drawText(Form("p_{T}>%d GeV", ptThr),xpos,ypos);
    drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy); 
    jumSun(-3.0,1,3.0,1);

    c1->SaveAs(Form("png/ES_eta_%s_%s_ptThr%d.png",sample.Data(),treePath.Data(),ptThr));
    c2[0]->SaveAs(Form("png/ES_etaDist_%s_%s_ptThr%d_%s.png",sample.Data(),treePath.Data(),ptThr,centSt[0]));
    c2[1]->SaveAs(Form("png/ES_etaDist_%s_%s_ptThr%d_%s.png",sample.Data(),treePath.Data(),ptThr,centSt[1]));
}
