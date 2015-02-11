#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <vector>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "TStyle.h"
#include "TSystem.h"

#include <TGraphAsymmErrors.h>
#include <TGraph.h>

#include <TMath.h>
#include <math.h>
#include "TFitResult.h"

#include <sstream>
#include <string>

#include "KYOcommonOpt.h"

using namespace std;
void formRapArr(Double_t binmin, Double_t binmax, string* arr);
void formPtArr(Double_t binmin, Double_t binmax, string* arr);

double fitERF(double *x, double *par);
double fitLine(double *x, double *par);

int draw_erfFit(char *stringA = "20140210_pt8bin", bool isPrompt=1)
{
	int contipar=2;

	//Read the single muon efficiency file from Kisoo
	//int nHist = 9;
	int nHist = 3;
	TFile* f_etabin[nHist];
	TH1D* h_etabin[nHist];
	for (Int_t i=0; i< nHist; i++){
		f_etabin[i] = new TFile(Form("./Closure_efficiency_total_etabin%dCS_1st_12bin_20140327.root",i+1));
		h_etabin[i] = (TH1D*)f_etabin[i]->Get("h3multi");
		cout << i <<"th h_etabin[i]" <<h_etabin[i] << endl;
	}

	gROOT->Macro("./JpsiStyle.C");
	gStyle->SetOptStat(0);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleSize(0.048,"xyz");

	string sampleType;

//	string etatitle[] = {"m240_m200", "m200_m150", "m150_m100", "m100_m050", "m050_m000", "m000_p050", "p050_p100", "p100_p150", "p150_p193"};
	string etatitle[] = {"m240_m080", "m080_p080", "p080_p193"};

	Double_t lowerLimit[3] = {0.8,3.3};

	// Set the erf(error function) formula to fit
	TF1* erfTmp[nHist]; 
	TF1* lineTmp[nHist]; 
	TF1* errFunc[nHist]; 
	//Double_t par[nHist][3];
	Double_t par[nHist][4]; // par0,1,2,3
	Double_t maxvalue[nHist]; // plateau for the highest bins 
	Double_t minvalue[nHist]; // plateau for the highest bins 

	//////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////
	///////////////////////// Fit AND Draw ///////////////////////////

	//latex box for the isPrompt, rapidity, pT info
	TLatex* latex = new TLatex();
	latex->SetNDC();
	latex->SetTextAlign(12);
	latex->SetTextSize(0.04);

	// --- Draw TGraphAsymmErrors
	TCanvas* c1 = new TCanvas("c1","c1",700,600);
  gStyle->SetOptTitle(0);
	//c1->SetLeftMargin(0.14);

	for (Int_t iy=0; iy<nHist; iy++) {
//		if (iy!=contipar) continue;
		minvalue[iy] = h_etabin[iy]->GetBinContent(1);
		maxvalue[iy] = h_etabin[iy]->GetBinContent(12);
		cout << "minvalue for " << etatitle[iy].c_str() << " = " << minvalue[iy] << endl;
		cout << "maxvalue for " << etatitle[iy].c_str() << " = " << maxvalue[iy] << endl;
		//erfTmp[iy] = new TF1(Form("erfTmp_%d",iy),fitERF,-40.0,40.0, 4);
		//errFunc[iy] = erfTmp[iy];
		errFunc[iy] = new TF1(Form("errFunc_%d",iy),fitERF,0.0,30.0,3);//0) plateau, 1) x-intrecept, 2) x broad 
		
		c1->cd();

		SetHistStyle(h_etabin[iy],5,0);
		h_etabin[iy]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_etabin[iy]->GetYaxis()->SetTitle("efficiency");
		h_etabin[iy]->GetXaxis()->SetRangeUser(0.,30.);
		h_etabin[iy]->GetYaxis()->SetRangeUser(0.,1.);
		h_etabin[iy]->GetXaxis()->CenterTitle();
		h_etabin[iy]->GetYaxis()->CenterTitle();
		h_etabin[iy]->Draw("ep");

		cout << "    "<< endl; 	
		cout << " ***** Fitting for :: "<<iy<< "th "<< h_etabin[iy]->GetName()<<endl;	
		cout << "maxvalue for " << etatitle[iy].c_str() << " = " << maxvalue[iy] << endl;
		cout << "maxvalue/2. = " << maxvalue[iy]/2. << endl;
		int counter =0;
		if (iy==0){
			errFunc[iy] = new TF1(Form("errFunc_%d",iy),fitERF,0.0,30.0,3);//0) plateau, 1) x-intrecept, 2) x broad 
			errFunc[iy]-> SetParameter(0,maxvalue[iy]);
			errFunc[iy]-> SetParLimits(0,maxvalue[iy]-0.005, maxvalue[iy]+0.002);
			errFunc[iy]-> SetParLimits(1,-0.7,0.7);
			errFunc[iy]-> SetParLimits(2,0.,3.);
		}
		else if (iy==1){
			h_etabin[iy]->SetBinContent(1,0);
			h_etabin[iy]->SetBinContent(2,0);
			errFunc[iy] = new TF1(Form("errFunc_%d",iy),fitERF,2.0,30.0,3);//0) plateau, 1) x-intrecept, 2) x broad 
			errFunc[iy]-> SetParameter(0,maxvalue[iy]);
			errFunc[iy]-> SetParLimits(0,maxvalue[iy]-0.02, maxvalue[iy]+0.001);
			errFunc[iy]-> SetParLimits(1,2,4);
			errFunc[iy]-> SetParLimits(2,2.,4.);
		}
		else if (iy==2){
			errFunc[iy] = new TF1(Form("errFunc_%d",iy),fitERF,0.0,30.0,3);//0) plateau, 1) x-intrecept, 2) x broad 
			errFunc[iy]-> SetParameter(0,maxvalue[iy]);
			errFunc[iy]-> SetParLimits(0,maxvalue[iy]-0.01, maxvalue[iy]+0.001);
			errFunc[iy]-> SetParLimits(1,-0.7,0.7);
			errFunc[iy]-> SetParLimits(2,1.5,3.);
		}

		TFitResultPtr res = h_etabin[iy]->Fit(Form("errFunc_%d",iy),"RSI");
		gPad->Update();
		if (0 != res->Status()) {
			cout << counter << " :-:-:KYO:-:-:Fitting is not converged :: " << h_etabin[iy]->GetName()<<endl;
			while (1) {
				counter++;
				res = h_etabin[iy]->Fit(Form("errFunc_%d",iy),"RSI");
				gPad->Update();
				if (0 != res->Status() || (counter >=4  && counter<8)) {
					errFunc[iy]-> SetParameter(0,maxvalue[iy]);
					errFunc[iy]->Update();
					gPad->Update();
				}
				if (0 == res->Status() || counter >= 8) { cout << "BREACK!!!!!!!!!!!!!!!"<<endl; break;}
			}
		}
		errFunc[iy]->GetParameters(&par[iy][0]);
		
		latex->DrawLatex(0.42, 0.35, Form("Status : %d (%d)",res->Status(),iy));
		latex->DrawLatex(0.42, 0.30, Form("par0 : %.3f",par[iy][0]));
		latex->DrawLatex(0.42, 0.26, Form("par1 : %.3f",par[iy][1]));
		latex->DrawLatex(0.42, 0.22, Form("par2 : %.3f",par[iy][2]));
		latex->DrawLatex(0.42, 0.17, etatitle[iy].c_str());
		//latex->DrawLatex(0.42, 0.16, Form("par3 : %.3f",par[iy][3]));
//		latex->DrawLatex(0.65, 0.30, sampleType.c_str());
//		latex->DrawLatex(0.65, 0.23, rapArr[iy].c_str());
//		latex->DrawLatex(0.65, 0.16, ptArr[iy].c_str());
		if (isPrompt==1) c1->SaveAs(Form("h_etabin_PR_%s.png",etatitle[iy].c_str()));
		else if (isPrompt==0) c1->SaveAs(Form("h_etabin_NP_%s.png",etatitle[iy].c_str()));
		c1->Clear();


	}

	// check the lowest values
	Double_t tmpx = 0;
	Double_t tmpy = 0;
	Double_t tmpx2 = 0;
	Double_t tmpy2 = 0;
	for (Int_t i=0; i<nHist; i++) {
//		if (i!=contipar) continue;
		if (i==0 || i==2) tmpx = lowerLimit[0];
		else tmpx = lowerLimit[1];
		tmpy = errFunc[i]->Eval(tmpx);
		//if (i==0) tmpx2=0.5;
		//tmpy2 = errFunc[i]->Eval(tmpx2);
		//if (tmpy <=0) {
			cout << i<<"th tmpx = " << tmpx << endl;
			cout << i<<"th tmpy = " << tmpy << endl;
		//}
	}

	// save as root file
	TFile *outFile = new TFile(Form("Eff_sigle_erfFit_%s_%s.root",sampleType.c_str(),stringA),"RECREATE");
	std::cout << "sampleType: " << sampleType.c_str() << std::endl;
	outFile->cd();
	for (Int_t i=0; i<nHist; i++) {
//		if (i!=contipar) continue;
		errFunc[i]->Write();
		h_etabin[i]->SetName(Form("h_etabin_%d",i));
		h_etabin[i]->Write();
	}
	outFile->Close();
	
	return 0;

}


void formRapArr(Double_t binmin, Double_t binmax, string* arr) {
	Double_t intMin, intMax; 
	Double_t fracMin = modf(binmin, &intMin);
	Double_t fracMax = modf(binmax, &intMax);
	if ( fracMin == 0 && fracMax == 0 ) {
		*arr = Form("%.0f<y_{lab}<%.0f", binmin, binmax);
	} else if ( fracMin != 0 && fracMax == 0 ) {
		*arr = Form("%.2f<y_{lab}<%.0f", binmin, binmax);
	} else if ( fracMin == 0 && fracMax != 0 ) {
		*arr = Form("%.0f<y_{lab}<%.2f", binmin, binmax);
	} else {
		*arr = Form("%.2f<y_{lab}<%.2f", binmin, binmax);
	}
}

void formPtArr(Double_t binmin, Double_t binmax, string* arr) {
	Double_t intMin, intMax; 
	Double_t fracMin = modf(binmin, &intMin);
	Double_t fracMax = modf(binmax, &intMax);
	if ( fracMin == 0 && fracMax == 0 ) {
		*arr = Form("%.0f<p_{T}<%.0f GeV/c", binmin, binmax);
	} else if ( fracMin != 0 && fracMax == 0 ) {
		*arr = Form("%.1f<p_{T}<%.0f GeV/c", binmin, binmax);
	} else if ( fracMin == 0 && fracMax != 0 ) {
		*arr = Form("%.0f<p_{T}<%.1f GeV/c", binmin, binmax);
	} else {
		*arr = Form("%.1f<p_{T}<%.1f GeV/c", binmin, binmax);
	}
}

double fitERF(double *x, double *par) {
	//return par[0]*TMath::Erf((x[0]-par[1])/par[2])+par[3];
	return par[0]*TMath::Erf((x[0]-par[1])/par[2]);
}

double fitLine(double *x, double *par) {
	return par[0]*x[0]+par[1];
}
