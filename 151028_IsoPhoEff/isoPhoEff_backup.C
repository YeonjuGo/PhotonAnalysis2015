/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */
#include <vector>
#include "TH1D.h"
#include "TStyle.h"
#include "TCut.h"

#include "../gedPhotonUtility.h"
static const long MAXTREESIZE = 10000000000;

const int nEtaCut_here = 4;
const float eta_i[nEtaCut_here] = {0.0, 0.0, 1.566, 2.0};
const float eta_f[nEtaCut_here] = {3.0, 1.4442, 2, 5.0};
const int nCentBin = 2;
const int cent_i[nCentBin] = {0,100};
const int cent_f[nCentBin] = {40,200};

void isoPhoEff(TString sample = "AllQCDPhoton30", TString treePath="ggHiNtuplizer",TString isoType = "sumIsoR",float isoCut=1.0, int rad=4, int evpOrder=2)
{
    const int nEP = 6;
    const double epmax = TMath::Pi();
    double epBin[nEP+1]; //{-epmax,-epmax/((double)nEP/2.),0.0,epmax/2.,epmax};
    for(int i=0; i<nEP+1;i++){
        double binWidth = (2*epmax)/nEP;
        epBin[i] = -epmax + binWidth*i;
    }
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    TString inputFileName = Form("../150908_voronoiStudy/skimFiles/jskim_%s_%s_inclusivePho_genIsoDR3_5GeV.root",sample.Data(),treePath.Data());
    TFile* fin = new TFile(inputFileName);
    TTree* tr = (TTree*)fin->Get("t_pho");

    TH1D* h1D_pt_den[5][5];//[eta][cent]
    TH1D* h1D_pt_num[5][5];//[eta][cent]
    TH1D* h1D_pt_eff[5][5];//[eta][cent]
    TH1D* h1D_dphi_den[5][5];//[eta][cent]
    TH1D* h1D_dphi_num[5][5];//[eta][cent]
    TH1D* h1D_dphi_eff[5][5];//[eta][cent]

    for(int ieta=0; ieta<nEtaCut_here; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            h1D_pt_den[ieta][icent] = new TH1D(Form("h1D_pt_den_eta%d_cent%d",ieta,icent),Form(";gen p_{T}^{#gamma};Isolation Efficiency(%s%d<%dGeV)",isoType.Data(),rad,(int)isoCut),10,0,100);
            h1D_pt_num[ieta][icent] = (TH1D*) h1D_pt_den[ieta][icent]->Clone(Form("h1D_pt_num_eta%d_cent%d",ieta,icent));
            h1D_pt_eff[ieta][icent] = (TH1D*) h1D_pt_den[ieta][icent]->Clone(Form("h1D_pt_eff_eta%d_cent%d",ieta,icent));
            h1D_dphi_den[ieta][icent] = new TH1D(Form("h1D_dphi_den_eta%d_cent%d",ieta,icent),Form(";#Delta#phi from the %d order Event Plane;Isolation Efficiency(%s%d<%dGeV)",evpOrder,isoType.Data(),rad,(int)isoCut),12,0,TMath::Pi());
            h1D_dphi_num[ieta][icent] = (TH1D*) h1D_dphi_den[ieta][icent]->Clone(Form("h1D_dphi_num_eta%d_cent%d",ieta,icent));
            h1D_dphi_eff[ieta][icent] = (TH1D*) h1D_dphi_den[ieta][icent]->Clone(Form("h1D_dphi_eff_eta%d_cent%d",ieta,icent));
        }
    }
    /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // Get Histograms
    TCanvas* ctmp = new TCanvas("ctmp", "", 100,100);
    for(int ieta=0; ieta<nEtaCut_here; ieta++){
        TString etaCut = Form("(phoEta>= %.3f && phoEta<%.3f)",eta_i[ieta],eta_f[ieta]);
        for(int icent=0;icent<nCentBin;icent++){
            TString centCut = Form("(hiBin>=%d && hiBin<%d)",cent_i[icent],cent_f[icent]);
            tr->Draw(Form("mcPt >>+ %s", h1D_pt_den[ieta][icent]->GetName()), Form("mcCalIsoDR0%d<5.0 && %s && %s",rad,etaCut.Data(),centCut.Data()));
            h1D_pt_den[ieta][icent] = (TH1D*)gDirectory->Get(h1D_pt_den[ieta][icent]->GetName());
            tr->Draw(Form("mcPt >>+ %s", h1D_pt_num[ieta][icent]->GetName()), Form("mcCalIsoDR0%d<5.0 && %s%d<%.1f && %s && %s",rad,isoType.Data(),rad,isoCut,etaCut.Data(),centCut.Data()));
            h1D_pt_num[ieta][icent] = (TH1D*)gDirectory->Get(h1D_pt_num[ieta][icent]->GetName());

            h1D_pt_eff[ieta][icent]->Divide(h1D_pt_num[ieta][icent],h1D_pt_den[ieta][icent],1,1,"B");
            tr->Draw(Form("phoDPhi_evtpl%d >>+ %s",evpOrder,h1D_dphi_den[ieta][icent]->GetName()), Form("mcCalIsoDR0%d<5 && %s && %s",rad,etaCut.Data(),centCut.Data()));
            h1D_dphi_den[ieta][icent] = (TH1D*)gDirectory->Get(h1D_dphi_den[ieta][icent]->GetName());
            tr->Draw(Form("phoDPhi_evtpl%d >>+ %s",evpOrder,h1D_dphi_num[ieta][icent]->GetName()), Form("mcCalIsoDR0%d<5.0 && %s%d<%.1f && %s && %s",rad,isoType.Data(),rad,isoCut,etaCut.Data(),centCut.Data()));
            h1D_dphi_num[ieta][icent] = (TH1D*)gDirectory->Get(h1D_dphi_num[ieta][icent]->GetName());
            h1D_dphi_eff[ieta][icent]->Divide(h1D_dphi_num[ieta][icent],h1D_dphi_den[ieta][icent],1,1,"B");
        }
    }

    /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // Cosmetics Hist
    for(int ieta=0; ieta<nEtaCut_here; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            h1D_pt_eff[ieta][icent]->SetMarkerStyle(20);
            h1D_dphi_eff[ieta][icent]->SetMarkerStyle(20);
            h1D_pt_eff[ieta][icent]->SetMarkerColor(2-icent);
            h1D_dphi_eff[ieta][icent]->SetMarkerColor(2-icent);
            h1D_pt_eff[ieta][icent]->GetYaxis()->SetRangeUser(0.,1.);
            h1D_dphi_eff[ieta][icent]->GetYaxis()->SetRangeUser(0.,1.);
        }
    }

    /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // Draw in Canvas

    double xpos = 0.5;
    double ypos = 0.8;
    double dy = 0.06;
    TLegend* leg = new TLegend(0.3,0.2,0.95,0.5);
    legStyle(leg);
    leg->AddEntry(h1D_pt_eff[0][0],"0-20 % (central)");
    leg->AddEntry(h1D_pt_eff[0][1],"50-100 % (peripheral)");

    TCanvas* c1 = new TCanvas("c1","eta tot",500,1000);
    c1->Divide(1,2);
    c1->cd(1);
    h1D_pt_eff[0][0]->Draw("p");
    h1D_pt_eff[0][1]->Draw("p same");
    drawText(Form("%.1f<|#eta|<%.1f",eta_i[0],eta_f[0]),xpos,ypos+2*dy);
    leg->Draw("same");
    c1->cd(2);
    h1D_dphi_eff[0][0]->Draw("p");
    h1D_dphi_eff[0][1]->Draw("p same");
    drawText(Form("%s",treePath.Data()),xpos,ypos-5*dy);
 
    TCanvas* c2 = new TCanvas("c2", "", 1000,600);
    c2->Divide(3,2);
    for(int ieta=1; ieta<nEtaCut_here; ieta++){
        c2->cd(ieta);
        h1D_pt_eff[ieta][0]->Draw("p");
        h1D_pt_eff[ieta][1]->Draw("p same");
        drawText(Form("%.2f<|#eta|<%.2f",eta_i[ieta],eta_f[ieta]),xpos,ypos+2*dy);
        c2->cd(ieta+3);
        h1D_dphi_eff[ieta][0]->Draw("p");
        h1D_dphi_eff[ieta][1]->Draw("p same");
        if(ieta==1) leg->Draw("same");
        if(ieta==2) drawText(Form("%s",treePath.Data()),xpos,ypos+2*dy);
        //if(ieta==3) drawText(Form("%s",isoType.Data()),xpos,ypos+2*dy);
    }
//    makeMultiPanelCanvas(can2,3,2,0.0,0.0,0.2,0.15,0.02);
   c1->SaveAs(Form("png/isoPhoEff_%s_%s_%s_R%d_evpOrder%d_etaDiff.png",sample.Data(),treePath.Data(),isoType.Data(),rad,evpOrder));
   c2->SaveAs(Form("png/isoPhoEff_%s_%s_%s_R%d_evpOrder%d_etaTot.png",sample.Data(),treePath.Data(),isoType.Data(),rad,evpOrder));
}

int main(){
    isoPhoEff();
    return 0;
}
