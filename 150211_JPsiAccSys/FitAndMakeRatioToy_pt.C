// latest version : 2015/02/06
// Yeonju Go 
// The code is based on FakeInput_Bplus() code by Hyunchul.
// The code fit the ratio(DATA/MC) histograms by Songkyo.
// pt(integrated of rapidity)

#include <iostream>
#include <iomanip>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TRandom.h>
#include "TFitResult.h"
#include <ctime>
#include "../HiForestAnalysis/hiForest.h"
#include "../gammaJetAnalysis/CutAndBinCollection2012.h"

Double_t ptArrNum[] = {0.0, 3.0, 4.0, 5.0, 6.5, 7.5, 8.5, 10., 14., 30.};
const Int_t nPt = sizeof(ptArrNum)/sizeof(Double_t)-1;

Double_t rapArrNumFB[] = {1.93, 1.5, 0.9, 0., -0.9, -1.5, -1.93, -2.4, -2.87};// for pt dist.
//Double_t rapArrNumBF[] = {-2.87, -2.4, -1.93, -1.5, -0.9, 0., 0.9, 1.5, 1.93};// for rap dist.
const Int_t nRap = sizeof(rapArrNumFB)/sizeof(Double_t)-1;


void FitAndMakeRatioToy_pt(int nToy=1, bool isPrompt=true, bool isPbp=true){
    gRandom->SetSeed(time(0));
    //	TFile* fin = new TFile(Form("/u/user/goyeonju/2015/pPbJPsiAnalysis/2015/004_closure/DataMcRecoIntegRap_8rap9pt/data_mc_reco_isPrompt%d.root",(int)isPrompt));//in KNU 
    TFile* fin = new TFile(Form("/home/goyeonju/CMS/2015/pPbJPsiAnalysis/2015/004_closure/DataMcRecoIntegRap_8rap9pt/data_mc_reco_isPrompt%d.root",(int)isPrompt));//in KNU 

    TH1D* hDataReco_Pbp=(TH1D*)fin->Get("hDataReco_Pbp");
    TH1D* hMCReco_Pbp=(TH1D*)fin->Get("hMCReco_Pbp");
    TH1D* hRatio_Pbp=(TH1D*)fin->Get("hRatio_Pbp");

    TH1D* hDataReco_pPb=(TH1D*)fin->Get("hDataReco_pPb");
    TH1D* hMCReco_pPb=(TH1D*)fin->Get("hMCReco_pPb");
    TH1D* hRatio_pPb=(TH1D*)fin->Get("hRatio_pPb");

    TFile* fout= new TFile(Form("toyResults_pt_number%d_isPrompt%d_isPbp%d.root",nToy,(int)isPrompt,(int)isPbp),"recreate");
    ////////////////////////////////////////////////////////////
    // Toy Tree setting
    double a1, a2, a3, a4;

    TTree* ditTree;
    ditTree = new TTree("ditTree","Toy MC fit parameters are stored in this tree");
    ditTree->Branch("a1",&a1,"a1/D");
    ditTree->Branch("a2",&a2,"a2/D");
    ditTree->Branch("a3",&a3,"a3/D");
    ditTree->Branch("a4",&a4,"a4/D");

    ////////////////////////////////////////////////////////////

    double Mean[nPt+1];
    double Err[nPt+1];
    double RatToy[nPt+1];
    for(int ipt=0;ipt<nPt;ipt++) {
        if(isPbp==true){
            Mean[ipt]=hRatio_Pbp->GetBinContent(ipt+1);
            Err[ipt]=hRatio_Pbp->GetBinError(ipt+1);
        } else{
            Mean[ipt]=hRatio_pPb->GetBinContent(ipt+1);
            Err[ipt]=hRatio_pPb->GetBinError(ipt+1);
        }	
    //    cout << "Bin " << ipt << " - " << Mean[ipt] << " + error " << Err[ipt] << endl;
    }
    TCanvas* c1 = new TCanvas("c1","",1200,1200);
    makeMultiPanelCanvas(c1,5,5,0.0,0.0,0.2,0.15,0.02);

    TH1D* htest[nPt];
    for(int ipt;ipt<nPt;ipt++){
        htest[ipt] = new TH1D(Form("htest_%d",ipt), Form("distribution of %dth toy(ratio) value",ipt), 50,0,2); 
    }

    for (int itoy = 0; itoy < 1000*nToy; itoy++){
        if (itoy% 2000 == 0)
            cout <<itoy<<" / "<<1000*nToy<<" "<<setprecision(2)<<(double)itoy/(1000*nToy)*100<<endl;	
        TH1D* hToyRatio = new TH1D(Form("hToyRatio%d",itoy), "Toy ratio", nPt,ptArrNum); 
        for (int ipt=0;ipt<nPt;ipt++) {
            RatToy[ipt]= gRandom->Gaus(Mean[ipt],Err[ipt]);
            //       printf("%.17f\n", RatToy[ipt]);

            hToyRatio -> SetBinContent(ipt+1,RatToy[ipt]);
            hToyRatio -> SetBinError(ipt+1,Err[ipt]);
            htest[ipt]->Fill(RatToy[ipt]);
        }
        TF1 *fRfitft;
        if(isPrompt){
            fRfitft= new TF1("fRfitft","[0]/(1+exp([1]*x))+[2]/x",0.0,30.0);//Prompt
            if(isPbp) fRfitft->SetParameters(8.98378e-1,-3.33945e-01,1.06557e+00);//Pbp ,isPbp=1
            else fRfitft->SetParameters(9.56389e-01,-2.95026e-01,9.19832e-01);//pPb, isPbp=0
            //fRfitft->SetParameters(0.95, -0.3, 0.92); 
            // fRfitft->SetParLimits(0,0.,2.);
            // fRfitft->SetParLimits(1,-2.,0.);
            // fRfitft->SetParLimits(2,0.90,0.94);//last curve
        } else {
            fRfitft= new TF1("fRfitft","[0]*TMath::Erf((x-[1])/[2])+[3]",0.0,30.0);//Non-prompt
            if(isPbp) fRfitft->SetParameters(1.86329e-01,5.64398e+00,2.76499e+00,9.38272e-01); //Pbp, isPbp=1
            else fRfitft->SetParameters(2.65086e-01,6.06224e+00,2.40510e+00,9.49416e-01); //pPb, isPbp=0
            //fRfitft->SetParameters(0.27, 6.06, 4.0, 0.95); //0) plateau, 1) x-shift, 2) x broad 3) y-shift
            fRfitft->SetParLimits(0,0.,0.7);
            fRfitft->SetParLimits(1,0.,10.);
            fRfitft->SetParLimits(2,0.,10.0);
            fRfitft->SetParLimits(3,0.,1.1);
        }
        hToyRatio->Fit("fRfitft","Q");

        if(itoy<25){
            c1->cd(itoy+1);
            hToyRatio->DrawCopy();
            fRfitft->DrawCopy();
            cout << hToyRatio->GetBinContent(1) << " , " << hToyRatio->GetBinContent(2)<< endl;
        }
        a1=fRfitft->GetParameter(0);
        a2=fRfitft->GetParameter(1);
        a3=fRfitft->GetParameter(2);
        if(isPrompt) a4=0.0;
        else a4=fRfitft->GetParameter(3);
        //     cout << "a1,a2,a3,a4 : " << a1 << ", " << a2 << ", " << a3 << ", " << a4 << endl;
        ditTree->Fill();

    } //itoy loop 

    fout->cd();	
    if(isPbp==true){
        hDataReco_Pbp->Write();
        hMCReco_Pbp->Write();
        hRatio_Pbp->Write();
    } else {
        hDataReco_pPb->Write();
        hMCReco_pPb->Write();
        hRatio_pPb->Write();
    }	
    //   fRfitft->Write();
    ditTree->Write();	
    for(int ipt;ipt<nPt;ipt++){
        htest[ipt]->Write();	
    }
    c1->Write();
    fout->Close();

}//end of main func.

