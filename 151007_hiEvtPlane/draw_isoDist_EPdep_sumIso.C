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

const int nCentBin = 2;
const int centBin[nCentBin+1] = {0,60,200};
const int nEP = 4;
const float epBin[nEP+1] = {-0.8,-0.4,0,0.4,0.8};

void draw_isoDist_EPdep_sumIso(TString outName="AllQCDPhoton30", bool isPromptPho=1, TString treePath="ggHiNtuplizer",int rad=2)
{
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    TString inputFileName = Form("../150908_voronoiStudy/skimFiles/jskim_%s_%s_isPromptPho%d.root",outName.Data(),treePath.Data(),(int)isPromptPho);
    TFile* fin = new TFile(inputFileName);
    //    fin->cd();
    TTree* tr = (TTree*)fin->Get("t_pho");    

    TH1D* h1D[5][2][4][4];//[eta][cent][evtpl][sumIsoType];
    int nbin = 100;
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<nEP;iep++){
                for(int isum=0;isum<3;isum++){
                    h1D[ieta][icent][iep][isum] = new TH1D(Form("h1D_eta%d_cent%d_ep%d_sumType%d",ieta,icent,iep,isum),";;",200,-100,500);
                }
            }
        }
    }

    TCanvas* ctmp=new TCanvas("ctmp","",300,300);
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<nEP;iep++){
                TCut etaCut = Form("(abs(phoEta)>= %.4f && abs(phoEta)< %.4f)", eta_gt[ieta], eta_lt[ieta]);
                TCut centCut = Form("(hiBin>= %d && hiBin< %d)", centBin[icent], centBin[icent+1]);
                TCut epCut = Form("(hiEvtPlane>= %.1f && hiEvtPlane< %.1f)", epBin[iep], epBin[iep+1]);
                TCut totalCut = etaCut && centCut && epCut;
                ctmp->cd();
                tr->Draw(Form("sumIsoR%d>>+%s",rad,h1D[ieta][icent][iep][0]->GetName()),totalCut);
                h1D[ieta][icent][iep][0]=(TH1D*)gDirectory->Get(h1D[ieta][icent][iep][0]->GetName()); 
                tr->Draw(Form("pfSumIso%d>>+%s",rad,h1D[ieta][icent][iep][1]->GetName()),totalCut);
                h1D[ieta][icent][iep][1]=(TH1D*)gDirectory->Get(h1D[ieta][icent][iep][1]->GetName()); 
                tr->Draw(Form("pfSumVsIso%d>>+%s",rad,h1D[ieta][icent][iep][2]->GetName()),totalCut);
                h1D[ieta][icent][iep][2]=(TH1D*)gDirectory->Get(h1D[ieta][icent][iep][2]->GetName()); 
            }
        }
    }

   double xpos = 0.5;
    double ypos = 0.8;
    double dy = 0.06;
    /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // AVG, RMS

    TLegend* leg2 = new TLegend(0.6,0.03,0.95,0.4);
    legStyle(leg2);
    TCanvas* can2 = new TCanvas("can2","", 1000,500);
    makeMultiPanelCanvas(can2,3,2,0.0,0.0,0.2,0.15,0.02);

    TH1D* h1D_avg[5][2][4];//[eta][cent][sumIsoType]
    TH1D* h1D_rms[5][2][4];//[eta][cent][sumIsoType]

    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int isum=0;isum<3;isum++){
                h1D_avg[ieta][icent][isum] = new TH1D(Form("h1D_avg_eta%d_cent%d_sumType%d",ieta,icent,isum),";hiEvtPlane;Avg",4,-TMath::Pi()/4.,TMath::Pi()/4.);
                h1D_rms[ieta][icent][isum] = new TH1D(Form("h1D_rms_eta%d_cent%d_sumType%d",ieta,icent,isum),";hiEvtPlane;RMS",4,-TMath::Pi()/4.,TMath::Pi()/4.);
                h1D_avg[ieta][icent][isum]->GetYaxis()->SetRangeUser(-20.0+rad*2.0,15.0+rad*30.0);
                h1D_rms[ieta][icent][isum]->GetYaxis()->SetRangeUser(-2.0+(rad*1.0),30.0+(rad*10.0));
                //h1D_rms[ieta][icent][isum]->GetYaxis()->SetRangeUser(0.0-((int)isPromptPho*5.0)+(rad*2.0),30.0-((int)isPromptPho*5.0)+(rad*2.0));
                h1D_avg[ieta][icent][isum]->SetMarkerStyle(20+isum);
                h1D_rms[ieta][icent][isum]->SetMarkerStyle(20+isum);
                h1D_avg[ieta][icent][isum]->SetMarkerSize(0.9);
                h1D_rms[ieta][icent][isum]->SetMarkerSize(0.9);
                h1D_avg[ieta][icent][isum]->SetMarkerColor(2*(icent+1));
                h1D_rms[ieta][icent][isum]->SetMarkerColor(2*(icent+1));
 
            }
/*          
            h1D_avg[ieta][icent][0]->SetMarkerColor(kRed+2*icent);
            h1D_rms[ieta][icent][0]->SetMarkerColor(kRed+2*icent);
            h1D_avg[ieta][icent][1]->SetMarkerColor(kGreen+2*icent);
            h1D_rms[ieta][icent][1]->SetMarkerColor(kGreen+2*icent);
            h1D_avg[ieta][icent][2]->SetMarkerColor(kBlue+2*icent);
            h1D_rms[ieta][icent][2]->SetMarkerColor(kBlue+2*icent); 
*/
        }
    }

    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            TString centSt;
            if(icent==0) centSt = "0-30%";
            else centSt = "30-100%";
            for(int iep=0;iep<nEP;iep++){
                for(int isum=0;isum<3;isum++){
                    h1D_avg[ieta][icent][isum]->SetBinContent(iep+1,h1D[ieta][icent][iep][isum]->GetMean());
                    h1D_avg[ieta][icent][isum]->SetBinError(iep+1,h1D[ieta][icent][iep][isum]->GetMeanError());
                    h1D_rms[ieta][icent][isum]->SetBinContent(iep+1,h1D[ieta][icent][iep][isum]->GetRMS());
                    h1D_rms[ieta][icent][isum]->SetBinError(iep+1,h1D[ieta][icent][iep][isum]->GetRMSError());
                }
            }
            if(ieta==1) {
                leg2->AddEntry(h1D_avg[ieta][icent][0], Form("oldSumIso%d %s",rad, centSt.Data() ));
                leg2->AddEntry(h1D_avg[ieta][icent][1], Form("pfSumIso%d %s",rad, centSt.Data() ));
                leg2->AddEntry(h1D_avg[ieta][icent][2], Form("pfVsSumIso%d %s",rad, centSt.Data() ));
            }

            for(int isum=0;isum<3;isum++){
                can2->cd(ieta);
                if(icent==0 && isum==0) h1D_avg[ieta][icent][isum]->Draw("pl");
                else h1D_avg[ieta][icent][isum]->Draw("pl same");
                can2->cd(ieta+3);
                if(icent==0 && isum==0) h1D_rms[ieta][icent][isum]->Draw("pl");
                else h1D_rms[ieta][icent][isum]->Draw("pl same");
            }
        }
        can2->cd(ieta);
        drawText(Form("%.2f<|eta|<%.2f",eta_gt[ieta],eta_lt[ieta]),xpos,ypos+2*dy); 
        if(ieta==1) leg2->Draw("same");
    } 
    can2->SaveAs(Form("png_sumIso/avg_rms_sumIsoTogether_R%d_%s_%s.png",rad,outName.Data(),treePath.Data()));       
}

int main(){
    draw_isoDist_EPdep_sumIso();
    return 0;
}
