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

Double_t ptArrNum[] = {0.0, 3.0, 4.0, 5.0, 6.5, 7.5, 8.5, 10., 14., 30.};
const Int_t nPt = sizeof(ptArrNum)/sizeof(Double_t)-1;

Double_t rapArrNumFB[] = {1.93, 1.5, 0.9, 0., -0.9, -1.5, -1.93, -2.4, -2.87};// for pt dist.
//Double_t rapArrNumBF[] = {-2.87, -2.4, -1.93, -1.5, -0.9, 0., 0.9, 1.5, 1.93};// for rap dist.
const Int_t nRap = sizeof(rapArrNumFB)/sizeof(Double_t)-1;


void FitAndMakeRatioToy_pt_test(int nToy=1, bool isPrompt=true, bool isPbp=true){
    gRandom->SetSeed(time(0));
    //	TFile* fin = new TFile(Form("/u/user/goyeonju/2015/pPbJPsiAnalysis/2015/004_closure/DataMcRecoIntegRap_8rap9pt/data_mc_reco_isPrompt%d.root",(int)isPrompt));//in KNU 
    TFile* fin = new TFile(Form("/home/goyeonju/CMS/2015/pPbJPsiAnalysis/2015/004_closure/DataMcRecoIntegRap_8rap9pt/data_mc_reco_isPrompt%d.root",(int)isPrompt));//in KNU 

    TH1D* hDataReco_Pbp=(TH1D*)fin->Get("hDataReco_Pbp");
    TH1D* hMCReco_Pbp=(TH1D*)fin->Get("hMCReco_Pbp");
    TH1D* hRatio_Pbp=(TH1D*)fin->Get("hRatio_Pbp");

    TH1D* hDataReco_pPb=(TH1D*)fin->Get("hDataReco_pPb");
    TH1D* hMCReco_pPb=(TH1D*)fin->Get("hMCReco_pPb");
    TH1D* hRatio_pPb=(TH1D*)fin->Get("hRatio_pPb");

    ////////////////////////////////////////////////////////////
    // Toy Tree setting
    Double_t a1, a2, a3, a4;

    TFile* fout= new TFile(Form("toyResults_pt_number%d_isPrompt%d_isPbp%d.root",nToy,(int)isPrompt,(int)isPbp),"recreate");
    TTree* ditTree;
    ditTree = new TTree("ditTree","Toy MC fit parameters are stored in this tree");
    ditTree->Branch("a1",&a1);
    ditTree->Branch("a2",&a2);
    ditTree->Branch("a3",&a3);
    ditTree->Branch("a4",&a4);

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
        //cout << "Bin " << ipt << " - " << Mean[ipt] << " + error " << Err[ipt] << endl;
    }

    TH1D* hToyRatio = new TH1D("hToyRatio", "Toy ratio", nPt,ptArrNum); 
    TH1D* htest[nPt];
    for(int ipt;ipt<nPt;ipt++){
        htest[ipt] = new TH1D(Form("htest_%d",ipt), Form("distribution of %dth toy(ratio) value",ipt), 50,0,2); 
    }

    //for (int itoy = 0; itoy < 10; itoy++){
        TF1 *fRfitft;
        if(isPrompt){
            fRfitft= new TF1("fRfitft","[0]/(1+exp([1]*x))+[2]/x",0.0,30.0);//Prompt
        } else {
            fRfitft= new TF1("fRfitft","[0]*TMath::Erf((x-[1])/[2])+[3]",0.0,30.0);//Non-prompt
        }
    for (int itoy = 0; itoy < 1000*nToy; itoy++){
        if (itoy% 2000 == 0)
            //cout <<itoy<<" / "<<1000*nToy<<" "<<setprecision(2)<<(double)itoy/(1000*nToy)*100<<endl;	
        hToyRatio->Clear();
        for (int ipt=0;ipt<nPt;ipt++) {
            RatToy[ipt]= gRandom->Gaus(Mean[ipt],Err[ipt]);
          //  printf("%.17f\n", RatToy[ipt]);

            hToyRatio -> SetBinContent(ipt+1,RatToy[ipt]);
            hToyRatio -> SetBinError(ipt+1,Err[ipt]);
            htest[ipt]->Fill(RatToy[ipt]);
        }
        if(isPrompt){
            if(isPbp) fRfitft->SetParameters(8.98378e-1,-3.33945e-01,1.06557e+00);//Pbp ,isPbp=1
            else fRfitft->SetParameters(9.56389e-01,-2.95026e-01,9.19832e-01);//pPb, isPbp=0
            //fRfitft->SetParameters(0.95, -0.3, 0.92); 
            // fRfitft->SetParLimits(0,0.,2.);
            // fRfitft->SetParLimits(1,-2.,0.);
            // fRfitft->SetParLimits(2,0.90,0.94);//last curve
        } else {
            if(isPbp) fRfitft->SetParameters(1.86329e-01,5.64398e+00,2.76499e+00,9.38272e-01); //Pbp, isPbp=1
            else fRfitft->SetParameters(2.65086e-01,6.06224e+00,2.40510e+00,9.49416e-01); //pPb, isPbp=0
            //fRfitft->SetParameters(0.27, 6.06, 4.0, 0.95); //0) plateau, 1) x-shift, 2) x broad 3) y-shift
            // fRfitft->SetParLimits(0,0.,1.);
            // fRfitft->SetParLimits(1,0.,20.);
            // fRfitft->SetParLimits(2,0.,20.);
            // fRfitft->SetParLimits(3,0.,2.);
        }
        /*        if(isPrompt==1 && isPbp==1)
                  fRfitft->SetParameters(8.98378e-1,-3.33945e-01,1.06557e+00);//Pbp ,isPbp=1
                  else if(isPrompt==1 && isPbp==0) 
                  fRfitft->SetParameters(9.56389e-01,-2.95026e-01,9.19832e-01);//pPb, isPbp=0
                  else if(isPrompt==0 && isPbp==1)
                  fRfitft->SetParameters(1.86329e-01,5.64398e+00,2.76499e+00,9.38272e-01);//Pbp
                  else if(isPrompt==0 && isPbp==0) 
                  fRfitft->SetParameters(2.65086e-01,6.06224e+00,2.40510e+00,9.49416e-01);//pPb
                  */

        hToyRatio->Fit("fRfitft", "Q");
#if 0
        //TFitResultPtr res = hToyRatio->Fit("fRfitft");
        TFitResultPtr res = hToyRatio->Fit("fRfitft","S");
        gPad->Update();
        if (0 != res->Status()) {
            cout << counter << " :-:-:Yeonju:-:-:Fitting is not converged :: " << hToyRatio->GetName()<<endl;
            while (1) {
                counter++;
                res = hToyRatio->Fit("fRfitft","RSI");
                gPad->Update();
                if (0 != res->Status() || (counter >=4  && counter<8)) {
                    fRfitft->Update();
                    gPad->Update();
                }
                if (0 == res->Status() || counter >= 8) { cout << "BREACK!!!!!!!!!!!!!!!"<<endl; break;}
            }
        }
#endif
        a1= (Double_t) fRfitft->GetParameter(0);
        a2= (Double_t) fRfitft->GetParameter(1);
        a3= (Double_t) fRfitft->GetParameter(2);
        if(isPrompt) a4=0.0;
        else a4= (Double_t) fRfitft->GetParameter(3);
        if (a3 < 0.00001)
             cout << "a1,a2,a3,a4 : " << a1 << ", " << a2 << ", " << a3 << ", " << a4 << endl;
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
    fout->Close();

}//end of main func.

