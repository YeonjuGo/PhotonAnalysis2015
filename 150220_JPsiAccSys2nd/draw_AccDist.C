// Author : Yeonju
// to draw Acceptance distributions

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>

const double shiftvar = -0.47; // conversion constant y=0(collision)==y=-0.47(LAB frame)  

void draw_AccDist(bool isPrompt=true, bool isPbp=true){
	using namespace std;

	// Definition of bin
	// --- pt Bin
	Double_t ptBinsArr[] = {0.0, 3.0, 4.0, 5.0, 6.5, 7.5, 8.5, 10.0, 14.0, 30.0}; // 8rap9pt
	//Double_t ptBinsArr[] = {5.0, 6.5, 10.0, 30.0}; // 6rap3pt
	const Int_t nPtBins = sizeof(ptBinsArr)/sizeof(double)-1;
	cout << "nPtBins=" << nPtBins << endl;

	Double_t AccCentArr[] = {0.087, 0.089, 0.131, 0.196, 0.288, 0.370, 0.447, 0.539, 0.690}; // Acceptance central value

    // --- y Bin //set to 1st run (For 2nd run, will be automatically changed later)
    Double_t yBinsArr[] = {-2.4, -1.97, -1.37, -0.47, 0.43, 1.03, 1.46, 1.93, 2.4}; // 8rap9pt
    //Double_t yBinsArr[] = {-2.4, -1.97, -1.37, -0.47, 0.43, 1.03, 1.46}; // 6rap3pt
    const Int_t nYBins = sizeof(yBinsArr)/sizeof(double)-1;
    cout << "nYBins=" << nYBins << endl;

	// for 2nd run
	Double_t yBinsArr2nd[nYBins+1] = {};
	for (Int_t i=0; i<nYBins+1; i++) {
		 yBinsArr2nd[i] = -1*yBinsArr[nYBins-i];
		cout <<"yBinsArr["<<i<<"] = " <<yBinsArr[i]<<endl;
		cout <<"yBinsArr2nd["<<i<<"] = " <<yBinsArr2nd[i]<<endl;
	}
	const Int_t nYBins2nd = sizeof(yBinsArr2nd)/sizeof(double)-1;

    ////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
    // Get the files and define histogrames

    TFile* fin[2][2];
    fin[0][0]=new TFile("AccDist_isPrompt0_isPbp0.root");
    fin[0][1]=new TFile("AccDist_isPrompt0_isPbp1.root");
    fin[1][0]=new TFile("AccDist_isPrompt1_isPbp0.root");
    fin[1][1]=new TFile("AccDist_isPrompt1_isPbp1.root");

    TH1D* hAccCompBin[2][2][nPtBins];
    for(int iprom=0;iprom<2;iprom++){
        for(int iPbp=0;iPbp<2;iPbp++){
            for(int ipt=0;ipt<nPtBins;ipt++){
                hAccCompBin[iprom][iPbp][ipt]=(TH1D*)fin[iprom][iPbp]->Get(Form("hAccCompBin%d",ipt));
                //new TH1D(Form("hAccCompBin%d_%d_%d",iprom,iPbp,ipt),Form("Acc. Dist. of %dth ptbin",ipt),6000, AccCentArr[ipt]-0.03, AccCentArr[ipt]+0.03);
                //hAccCompBin[iprom][iPbp][ipt]->Sumw2();
            }
        }
    }

    TCanvas* c1[2][2];
    for(int iprom=0;iprom<2;iprom++){
        for(int iPbp=0;iPbp<2;iPbp++){
            c1[iprom][iPbp]=new TCanvas(Form("c1_%d_%d",iprom,iPbp),Form("AccDist_isPrompt%d_isPbp%d",iprom,iPbp), 1200,900);
            c1[iprom][iPbp]->Divide(4,3);
            for(int ipt=0;ipt<nPtBins;ipt++){
                c1[iprom][iPbp]->cd(ipt+1);
                hAccCompBin[iprom][iPbp][ipt]->Draw("hist"); 
                hAccCompBin[iprom][iPbp][ipt]->GetXaxis()->SetRangeUser(AccCentArr[ipt]-0.03,AccCentArr[ipt]+0.03); 
            }
        }
    }
    
    TFile* outf = new TFile("draw_AccDist.root","recreate");
    outf->cd();
    c1[0][0]->Write();
    c1[1][0]->Write();
    c1[0][1]->Write();
    c1[1][1]->Write();
    outf->Close();
}
