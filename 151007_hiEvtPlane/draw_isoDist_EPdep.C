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

void draw_isoDist_EPdep(TString outName="AllQCDPhoton30", bool isPromptPho=1, TString treePath="ggHiNtuplizer", TString var="pho_ecalClusterIsoR2", int evtplNum=2, double xmin = -100, double xmax = 400, int nbin = 100)
{
    int evpOrder = 1;
    if(evtplNum==2) evpOrder = 1;
    else if(evtplNum==8) evpOrder = 2;
    else if(evtplNum==15) evpOrder = 3;
    else if(evtplNum==21) evpOrder = 4;
    const int nEP = 4;
    const double epmax = TMath::Pi()/(double)evpOrder;
    double epBin[nEP+1]; //{-epmax,-epmax/((double)nEP/2.),0.0,epmax/2.,epmax};
    for(int i=0; i<nEP+1;i++){
        double binWidth = (2*epmax)/(double)nEP;
        epBin[i] = -epmax + binWidth*i;
    }
    TString promSt = "";
    if(isPromptPho) promSt = "Prmopt photon";
    else promSt = "Non-prompt photon";
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleSize(0.05,"X");
    TString inputFileName = Form("../150908_voronoiStudy/skimFiles/jskim_%s_%s_isPromptPho%d_hiEvtPlane%d.root",outName.Data(),treePath.Data(),(int)isPromptPho,evtplNum);
    TFile* fin = new TFile(inputFileName);
    TTree* tr = (TTree*)fin->Get("t_pho");    
    TCanvas* can = new TCanvas("can","", 1000,600);
    can->Divide(3,2);
    TH1D* h1D[5][5][10];//[eta][cent][evtpl]
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<nEP;iep++){
                h1D[ieta][icent][iep] = new TH1D(Form("h1D_eta%d_cent%d_ep%d",ieta,icent,iep),Form(";%s;",var.Data()),nbin,xmin,xmax);
            }
        }
    }

    TCanvas* ctmp=new TCanvas("ctmp","",300,300);
    cout << tr << endl;
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<nEP;iep++){
                TCut etaCut = Form("(abs(phoEta)>= %.4f && abs(phoEta)< %.4f)", eta_gt[ieta], eta_lt[ieta]);
                TCut centCut = Form("(hiBin>= %d && hiBin< %d)", centBin[icent], centBin[icent+1]);
                TCut epCut = Form("(hiEvtPlane>= %.1f && hiEvtPlane< %.1f)", epBin[iep], epBin[iep+1]);
                TCut totalCut = etaCut && centCut && epCut;
                ctmp->cd();
                //tr->Draw(Form("%s>>+h1D_eta%d_cent%d_ep%d",var.Data(),ieta,icent,iep),totalCut);
                tr->Draw(Form("%s>>+%s",var.Data(),h1D[ieta][icent][iep]->GetName()),totalCut);
                h1D[ieta][icent][iep]=(TH1D*)gDirectory->Get(h1D[ieta][icent][iep]->GetName()); 
            }
        }
    }

    TLegend* leg = new TLegend(0.4,0.2,0.95,0.6);
    legStyle(leg);
    double xpos = 0.5;
    double ypos = 0.8;
    double dy = 0.06;
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<nEP;iep++){
                can->cd(ieta+3*icent);
                gPad->SetLogy();
                h1D[ieta][icent][iep]->SetLineColor(2*iep+2);

                if(ieta==1 && icent==0) leg->AddEntry(h1D[ieta][icent][iep],Form("%.1f < EP[%d] < %.1f",epBin[iep],evpOrder,epBin[iep+1]));
                if(iep==0) h1D[ieta][icent][iep]->Draw("hist");
                else h1D[ieta][icent][iep]->Draw("same hist");
            }
            if(ieta==1){
                can->cd(ieta+3*icent);
                drawText(Form("%d < hiBin < %d", centBin[icent], centBin[icent+1]), xpos, ypos-dy);
            }
        }//cent loop

        if(ieta==1) { can->cd(ieta); leg->Draw("same"); }
        else if(ieta==3) { 
            can->cd(ieta);
            drawText(Form("%s",promSt.Data()),xpos,ypos-2*dy); 
            drawText(Form("%s",treePath.Data()),xpos,ypos-3*dy); 
        } 
        can->cd(ieta);
        drawText(Form("%.2f<|#eta|<%.2f",eta_gt[ieta],eta_lt[ieta]),xpos,ypos); 
    }//eta loop
    can->SaveAs(Form("png_logScale/hiEvtPlanes%d_%s_isPromptPho%d_%s.png",evtplNum, var.Data(),(int)isPromptPho,treePath.Data())); 
    //can->SaveAs(Form("png/%s_%s_%s.png",var.Data(),outName.Data(),treePath.Data()));       

    /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // AVG, RMS
#if 0
    TLegend* leg2 = new TLegend(0.5,0.1,0.95,0.5);
    legStyle(leg2);
    TCanvas* can2 = new TCanvas("can2","", 1000,500);
    makeMultiPanelCanvas(can2,3,2,0.0,0.0,0.2,0.15,0.02);

    TH1D* h1D_avg_temp = new TH1D("h1D_avg_temp","",4,-TMath::Pi()/2.,TMath::Pi()/2.);
    TH1D* h1D_avg[5][2];//[eta][cent][evtpl]
    TH1D* h1D_rms[5][2];//[eta][cent][evtpl]

    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            h1D_avg[ieta][icent] = new TH1D(Form("h1D_avg_eta%d_cent%d",ieta,icent),";hiEvtPlane;Avg",4,-TMath::Pi()/4.,TMath::Pi()/4.);
            h1D_rms[ieta][icent] = new TH1D(Form("h1D_rms_eta%d_cent%d",ieta,icent),";hiEvtPlane;RMS",4,-TMath::Pi()/4.,TMath::Pi()/4.);
            h1D_avg[ieta][icent]->GetYaxis()->SetRangeUser(-10.0,10.0);
            h1D_rms[ieta][icent]->GetYaxis()->SetRangeUser(0.0,10.0);
            h1D_avg[ieta][icent]->SetMarkerStyle(20);
            h1D_rms[ieta][icent]->SetMarkerStyle(20);
            h1D_avg[ieta][icent]->SetMarkerColor(2*(icent+1));
            h1D_rms[ieta][icent]->SetMarkerColor(2*(icent+1));
        }
    }
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<nEP;iep++){
                h1D_avg[ieta][icent]->SetBinContent(iep+1,h1D[ieta][icent][iep]->GetMean());
                h1D_avg[ieta][icent]->SetBinError(iep+1,h1D[ieta][icent][iep]->GetMeanError());
                h1D_rms[ieta][icent]->SetBinContent(iep+1,h1D[ieta][icent][iep]->GetRMS());
                h1D_rms[ieta][icent]->SetBinError(iep+1,h1D[ieta][icent][iep]->GetRMSError());
            }
        if(ieta==1) leg2->AddEntry(h1D_avg[ieta][icent], Form("%d < hiBin < %d", centBin[icent], centBin[icent+1]));
        can2->cd(ieta);
        if(icent==0) h1D_avg[ieta][icent]->Draw("pl");
        else h1D_avg[ieta][icent]->Draw("pl same");
        can2->cd(ieta+3);
        if(icent==0) h1D_rms[ieta][icent]->Draw("pl");
        else h1D_rms[ieta][icent]->Draw("pl same");
        }
    can2->cd(ieta);
    drawText(Form("%.2f<|eta|<%.2f",eta_gt[ieta],eta_lt[ieta]),xpos,ypos-dy); 
    drawText(Form("%s",var.Data()),xpos,ypos); 
    if(ieta==1) leg2->Draw("same");
    } 
    can2->SaveAs(Form("png_logScale/avg_rms_%s_%s_%s.png",var.Data(),outName.Data(),treePath.Data()));       
#endif
}

int main(){
    draw_isoDist_EPdep();
    return 0;
}
