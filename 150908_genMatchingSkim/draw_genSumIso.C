/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */
#include <vector>
#include "TH2D.h"
#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#include "../gedPhotonUtility.h" 
static const long MAXTREESIZE = 10000000000;

void draw_genSumIso(TString outName="EmEnrichedDijet30", TString treePath="ggHiNtuplizer", float ptThr=20, TString cond="abs(phoEta)<=1.44")
{
    TString inputFileName = Form("skimFiles/jskim_%s_%s.root",outName.Data(),treePath.Data());
    TFile* fin = new TFile(inputFileName);
    fin->cd();
    TTree* tr = (TTree*)fin->Get("t_phogen");
    TCanvas* can[5];
    TH2D* h2D[5][3];
    TH1D* h1D;
    for(int i=0; i<7; i++){
       for(int j=0;j<3;j++){
            h2D[i][j] = new TH2D(Form("h2D_rad%d_iso%d",i,j),";;genSumIso",35,-50,300,50,0,500);
        }
    }

    TCanvas* ctmp=new TCanvas("ctmp","",300,300);
    cout << tr << endl;
    for(int i=2; i<6; i++){
        ctmp->cd();
        tr->Draw(Form("genSumIso%d:pho_ecalClusterIsoR%d+pho_hcalRechitIsoR%d+pho_trackIsoR%dPtCut20>>+%s",i,i,i,i,h2D[i][0]->GetName()),cond.Data(),"colz");
        h2D[i][0]=(TH2D*)gDirectory->Get(h2D[i][0]->GetName()); 
        tr->Draw(Form("genSumIso%d:pfcIso%d+pfpIso%d+pfnIso%d>>+%s",i,i,i,i,h2D[i][1]->GetName()),cond,"colz");
        h2D[i][1]=(TH2D*)gDirectory->Get(h2D[i][1]->GetName()); 
        tr->Draw(Form("genSumIso%d:pfcVsIso%d+pfpVsIso%d+pfnVsIso%d>>+%s",i,i,i,i,h2D[i][2]->GetName()),cond,"colz");
        h2D[i][2]=(TH2D*)gDirectory->Get(h2D[i][2]->GetName()); 
    }
   
    for(int i=2; i<6; i++){ 
        can[i] = new TCanvas(Form("c%d",i),"", 1000,300);
        can[i]->Divide(3,1);
        can[i]->cd(1);
        h2D[i][0]->Draw("colz");
        can[i]->cd(2);
        h2D[i][1]->Draw("colz");
        can[i]->cd(3);
        h2D[i][2]->Draw("colz");
        can[i]->SaveAs(Form("png/2D_genSumIso_recoSumIso_%s_%s_%s_rad%d.png",outName.Data(),treePath.Data(),cond.Data(),i));       
    }    
   // fin->Close();
}
