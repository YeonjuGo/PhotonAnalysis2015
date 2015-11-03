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

void draw_isoDist_EPdep_inclusivePho(TString treePath="ggHiNtuplizer", TString var="pho_ecalClusterIsoR2", int evtplNum=2, double xmin = -50, double xmax = 150, int nbin = 50)
{
    int evpOrder = 1;
    if(evtplNum==2) evpOrder = 1;
    else if(evtplNum==8) evpOrder = 2;
    else if(evtplNum==15) evpOrder = 3;
    else if(evtplNum==21) evpOrder = 4;
    const int nEP = 4;
    const double epmax = TMath::Pi();
    double epBin[nEP+1]={0.0,epmax/nEP,2*epmax/nEP,3*epmax/nEP,epmax};
/*    for(int i=0; i<nEP+1;i++){
        double binWidth = (2*epmax)/(double)nEP;
        epBin[i] = -epmax + binWidth*i;
    }
*/
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleSize(0.05,"X");
    //gStyle->SetLineWidth(1.3);

    TString inputFileName = Form("../150908_voronoiStudy/skimFiles/jskim_%s_inclusivePho_genIsoDR3_5GeV.root",treePath.Data());
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
    //cout << tr << endl;
    for(int ieta=1; ieta<nEtaCut; ieta++){
        for(int icent=0;icent<nCentBin;icent++){
            for(int iep=0;iep<nEP;iep++){
                TCut etaCut = Form("(abs(phoEta)>= %.4f && abs(phoEta)< %.4f)", eta_gt[ieta], eta_lt[ieta]);
                TCut centCut = Form("(hiBin>= %d && hiBin< %d)", centBin[icent], centBin[icent+1]);
                TCut epCut = Form("(abs(phoDPhi_evtpl%d) >= %.1f  && abs(phoDPhi_evtpl%d) < %.1f)",evpOrder, epBin[iep], evpOrder, epBin[iep+1]);
                TCut hoeCut = "(phoHoverE<0.1)"; 
                TCut totalCut = etaCut && centCut && epCut && hoeCut;
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
            double yRange=0.0;
            for(int iep=0;iep<nEP;iep++){
                can->cd(ieta+3*icent);
//                gPad->SetLogy();
                h1D[ieta][icent][iep]->SetLineColor(2*iep+2);

                if(ieta==1 && icent==0) leg->AddEntry(h1D[ieta][icent][iep],Form("%.1f<#Delta#phi(EP[%d])<%.1f",epBin[iep],evpOrder,epBin[iep+1]));
                if(iep==0) h1D[ieta][icent][iep]->Draw("hist");
                else h1D[ieta][icent][iep]->Draw("same hist");
                double tmpRange = cleverRange(h1D[ieta][icent][iep]);
                if(yRange<tmpRange) yRange=tmpRange;
            }
            if(ieta==1){
                can->cd(ieta+3*icent);
                drawText(Form("%d < hiBin < %d", centBin[icent], centBin[icent+1]), xpos, ypos-dy);
            }
            
            for(int iep=0;iep<nEP;iep++){
               h1D[ieta][icent][iep]->SetAxisRange(1.e-3,yRange,"Y");
            }
        }//cent loop
        if(ieta==1) { can->cd(ieta); leg->Draw("same"); }
        else if(ieta==3) { 
            can->cd(ieta);
            drawText(Form("Inclusive Photon"),xpos,ypos-2*dy); 
            drawText(Form("%s",treePath.Data()),xpos,ypos-3*dy); 
        } 
        can->cd(ieta);
        drawText(Form("%.2f<|#eta|<%.2f",eta_gt[ieta],eta_lt[ieta]),xpos,ypos); 
    }//eta loop
    can->SaveAs(Form("png/EPorder%d_%s_inclusivePho_%s.png",evpOrder, var.Data(),treePath.Data())); 
    //can->SaveAs(Form("png_logScale/hiEvtPlanes%d_%s_inclusivePho_%s.png",evtplNum, var.Data(),treePath.Data())); 
}

int main(){
    draw_isoDist_EPdep_inclusivePho();
    return 0;
}
