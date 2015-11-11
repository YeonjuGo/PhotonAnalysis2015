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

void isoPhoEff_isoDep_noGenIsoCut(TString sample = "EmEnrichedDijet30", TString treePath="ggHiNtuplizer",TString isoType = "sumIsoR", int rad=4, int evpOrder=2)
{
    TString var = Form("%s%d",isoType.Data(),rad);

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
    TString inputFileName = Form("../150908_voronoiStudy/skimFiles/jskim_%s_%s_genMatched_inclusivePho.root",sample.Data(),treePath.Data());
    TFile* fin = new TFile(inputFileName);
    TTree* tr = (TTree*)fin->Get("t_pho");

    TH1D* h1D_iso_eff[5][5];//[eta][cent]

    const int nBin = 20;
    const double xmin = 0;
    const double xmax = 100;
    double binWidth = (xmax-xmin)/nBin;
    double middle_0thBin = binWidth/2.;
    for(int ieta=0; ieta<nEtaCut_here; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            //iso
            h1D_iso_eff[ieta][icent] = new TH1D(Form("h1D_iso_eff_eta%d_cent%d",ieta,icent),Form(";%s;Isolation Efficiency",var.Data()),nBin,xmin,xmax);
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
            for(int j=0;j<nBin;j++){
                double den = tr->GetEntries(Form("%s && %s", etaCut.Data(), centCut.Data()));
                double num = tr->GetEntries(Form("%s>-100 && %s< %.1f && %s && %s",var.Data(), var.Data(), middle_0thBin+j*binWidth, etaCut.Data(), centCut.Data()));
                double eff = num/den;
                h1D_iso_eff[ieta][icent]->SetBinContent(j+1,eff);
            }
        }
    }

    /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // Cosmetics Hist
    for(int ieta=0; ieta<nEtaCut_here; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            h1D_iso_eff[ieta][icent]->SetMarkerStyle(20);
            h1D_iso_eff[ieta][icent]->SetMarkerColor(2-icent);
            h1D_iso_eff[ieta][icent]->GetYaxis()->SetRangeUser(0.,1.);
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
    leg->AddEntry(h1D_iso_eff[0][0],"0-20 % (central)");
    leg->AddEntry(h1D_iso_eff[0][1],"50-100 % (peripheral)");

    TCanvas* c1 = new TCanvas("c1","eta tot",500,500);
    c1->cd();
    gPad->SetGrid();
    h1D_iso_eff[0][0]->Draw("p");
    h1D_iso_eff[0][1]->Draw("p same");
    drawText(Form("%.1f<|#eta|<%.1f",eta_i[0],eta_f[0]),xpos,ypos+2*dy);
    leg->Draw("same");
    drawText(Form("%s",treePath.Data()),xpos,ypos-5*dy);

    TCanvas* c2 = new TCanvas("c2", "", 1000,300);
    c2->Divide(3,1);
    for(int ieta=1; ieta<nEtaCut_here; ieta++){
        c2->cd(ieta);
        gPad->SetGrid();
        h1D_iso_eff[ieta][0]->Draw("p");
        h1D_iso_eff[ieta][1]->Draw("p same");
        drawText(Form("%.2f<|#eta|<%.2f",eta_i[ieta],eta_f[ieta]),xpos,ypos+2*dy);
        if(ieta==1) leg->Draw("same");
        if(ieta==2) drawText(Form("%s",treePath.Data()),xpos,ypos-5*dy);
        //if(ieta==3) drawText(Form("%s",isoType.Data()),xpos,ypos+2*dy);
    }
    c1->SaveAs(Form("png/isoPhoEff_isoDep_noGenIsoCut_%s_%s_%s_R%d_evpOrder%d_etaTot.png",sample.Data(),treePath.Data(),isoType.Data(),rad,evpOrder));
    c2->SaveAs(Form("png/isoPhoEff_isoDep_noGenIsoCut_%s_%s_%s_R%d_evpOrder%d_etaDiffrential.png",sample.Data(),treePath.Data(),isoType.Data(),rad,evpOrder));
}

int main(){
    isoPhoEff_isoDep_noGenIsoCut();
    return 0;
}
