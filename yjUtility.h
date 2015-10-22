#ifndef UTILITY_Yeonju_H
#define UTILITY_Yeonju_H
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TDatime.h>
#include <iostream>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TStyle.h>

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>
using namespace std;

const int MAXGENPARTICLES = 50000;  // number of gen particles can be large
const int MAXPHOTONS = 500;
const int PDG_PHOTON = 22;
const int cutmcStatus = 1;
const double cutdeltaR = 0.2; // 0.3    // cut for matching gen and reco. particles
//const float cutetaBarrel = 1.4791;                      // cut to separate photons into Barrel and Endcap photons
const float cutetaBarrel = 1.4442;                      // cut to separate photons into Barrel and Endcap photons
const float cutetaBarrelGap = 1.566;                      // cut to separate photons into Barrel and Endcap photons
const float cutetaEndCap = 2;                           // cut to separate photons in Endcap into 2.
const int cutmcMomPID_pi0 = 111;

const int nEtaCut = 4;
const int nCentCut = 5;
const int nMomIdCut = 5;
const int nRadius = 6;

float eta_gt[nEtaCut] = {  0.0,            0.0, cutetaBarrelGap, cutetaEndCap};
float eta_lt[nEtaCut] = {5.0, cutetaBarrel, cutetaEndCap,       5.0};
int hiBin_gt[nCentCut] = {-999,  0, 20,  60, 100};
int hiBin_lt[nCentCut] = { 999, 20, 60, 100, 200};
int mcMomPID_gt[nMomIdCut] = {-999, 21, -999,  22, 110};
int mcMomPID_lt[nMomIdCut] = { 999, 23,   22, 999, 112};

enum condition { noC, hoeC, sigmaC, hoeAndSigmaC };
TString getCondSuffix ( condition cond_) {
  if (cond_ == noC) return "";
  if (cond_ == hoeC) return "_hoe0.1";
  if (cond_ == sigmaC) return "_sigma0.01";
  if (cond_ == hoeAndSigmaC) return "_hoe0.1_sigma0.01";
  return "NULL";
}
TString getCondDirName ( condition cond_) {
  if (cond_ == noC) return "_noCut";
  if (cond_ == hoeC ) return "_hoeC";
  if (cond_ == sigmaC ) return "_sigmaC";
  if (cond_ == hoeAndSigmaC) return "_hoeAndSigmaC";
  return "NULL";
}

TString getCondSuffix ( int cond_) {
  if (cond_ == noC) return "";
  if (cond_ == hoeC) return "_hoe0.1";
  if (cond_ == sigmaC) return "_sigma0.01";
  if (cond_ == hoeAndSigmaC) return "_hoe0.1_sigma0.01";
  return "NULL";
}
TString getCondDirName ( int cond_) {
  if (cond_ == noC) return "_noCut";
  if (cond_ == hoeC ) return "_hoeC";
  if (cond_ == sigmaC ) return "_sigmaC";
  if (cond_ == hoeAndSigmaC) return "_hoeAndSigmaC";
  return "NULL";
}

void legStyle( TLegend *a=0 , TString head="")
{
  a->SetBorderSize(0);
  a->SetHeader(head);
//  a->SetTextFont(62);
//  a->SetTextSize(17);
//  a->SetLineColor(1);
//  a->SetLineStyle(1);
//  a->SetLineWidth(1);
//  a->SetFillColor(0);
  a->SetFillStyle(0);

}
void graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
{
	g1->SetLineStyle(lstyle);
	g1->SetLineColor(lcolor);
	g1->SetLineWidth(lwidth);
	g1->SetMarkerStyle(mstyle);
	g1->SetMarkerColor(mcolor);
	g1->SetMarkerSize(msize);
}

void hLineStyle(TH1 *h1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t lwidth=1, Int_t lfst=0, Int_t lfcolor=0)
{
	h1->SetLineStyle(lstyle);
	h1->SetLineColor(lcolor);
	h1->SetLineWidth(lwidth);
	h1->SetFillStyle(lfst);
	h1->SetFillColor(lfcolor);
}

void hMarkerStyle(TH1 *h1=0, Int_t mstyle=20, Int_t mcolor=1, Double_t msize=1.0)
{
	h1->SetMarkerStyle(mstyle);
	h1->SetMarkerColor(mcolor);
	h1->SetMarkerSize(msize);
}
void drawText(const char *text, float xp, float yp, int textColor=kBlack, int textSize=14){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(43);
	//   if(bold)tex->SetTextFont(43);
	tex->SetTextSize(textSize);
	tex->SetTextColor(textColor);
	tex->SetLineWidth(1);
	tex->SetNDC();
	tex->Draw();
}
void jumSun(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
	TLine* t1 = new TLine(x1,y1,x2,y2);
	t1->SetLineWidth(width);
	t1->SetLineStyle(7);
	t1->SetLineColor(color);
	t1->Draw();
}
double findCross(TH1* h1, TH1* h2, double& frac, double& effi, double& fracErr, double& effiErr){
	Int_t nBins = h1->GetNbinsX();
	double crossVal =0;
	int binAt0 = h1->FindBin(0);
	for(Int_t ix=binAt0; ix<=nBins ;ix++){
		float yy1 = h1->GetBinContent(ix);
		float yy2 = h2->GetBinContent(ix);
		if(yy2>yy1) {
			crossVal= h1->GetBinLowEdge(ix);
			break;
		}
	}
	int crossBin = h1->FindBin(crossVal);
	frac = 1 - (h2->Integral(1,crossBin) / h1->Integral(1,crossBin) );
	effi = ( h1->Integral(1,crossBin) / h1->Integral() );
	fracErr = frac * TMath::Sqrt( (1./h2->Integral(1,crossVal)) + (1./h1->Integral(1,crossVal)) );
	effiErr = ( TMath::Sqrt(h1->Integral(1,crossVal)) / h1->Integral() ) * TMath::Sqrt(1 - (h1->Integral(1,crossVal)/h1->Integral()) );

	return crossVal;
}

void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns,
		const Int_t rows, const Float_t leftOffset=0.,
		const Float_t bottomOffset=0.,
		const Float_t leftMargin=0.2,
		const Float_t bottomMargin=0.2,
		const Float_t edge=0.05) {
	if (canv==0) {
		//Error("makeMultiPanelCanvas","Got null canvas.");
		return;
	}
	canv->Clear();

	TPad* pad[columns][rows];

	Float_t Xlow[columns];
	Float_t Xup[columns];
	Float_t Ylow[rows];
	Float_t Yup[rows];
	Float_t PadWidth =
		(1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
				(1.0/(1.0-edge))+(Float_t)columns-2.0);
	Float_t PadHeight =
		(1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
				(1.0/(1.0-edge))+(Float_t)rows-2.0);
	Xlow[0] = leftOffset;
	Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
	Xup[columns-1] = 1;
	Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

	Yup[0] = 1;
	Ylow[0] = 1.0-PadHeight/(1.0-edge);
	Ylow[rows-1] = bottomOffset;
	Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

	for(Int_t i=1;i<columns-1;i++) {
		Xlow[i] = Xup[0] + (i-1)*PadWidth;
		Xup[i] = Xup[0] + (i)*PadWidth;
	}
	Int_t ct = 0;
	for(Int_t i=rows-2;i>0;i--) {
		Ylow[i] = Yup[rows-1] + ct*PadHeight;
		Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
		ct++;
	}
	TString padName;
	for(Int_t i=0;i<columns;i++) {
		for(Int_t j=0;j<rows;j++) {
			canv->cd();
			padName = Form("p_%d_%d",i,j);
			pad[i][j] = new TPad(padName.Data(),padName.Data(),
					Xlow[i],Ylow[j],Xup[i],Yup[j]);
			if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
			else pad[i][j]->SetLeftMargin(0);

			if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
			else pad[i][j]->SetRightMargin(0);

			if(j==0) pad[i][j]->SetTopMargin(edge);
			else pad[i][j]->SetTopMargin(0);

			if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
			else pad[i][j]->SetBottomMargin(0);

			pad[i][j]->Draw();
			pad[i][j]->cd();
			pad[i][j]->SetNumber(columns*j+i+1);
		}
	}
}



Double_t getDPHI( Double_t phi1, Double_t phi2) {
        Double_t dphi = phi1 - phi2;

        if ( dphi > 3.141592653589 )
                dphi = dphi - 2. * 3.141592653589;
        if ( dphi <= -3.141592653589 )
                dphi = dphi + 2. * 3.141592653589;

        if ( TMath::Abs(dphi) > 3.141592653589 ) {
                std::cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << std::endl;
        }

        return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2){
        return eta1 - eta2;
}

Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
        Double_t theDphi = getDPHI( phi1, phi2);
        Double_t theDeta = eta1 - eta2;
        return TMath::Sqrt ( theDphi*theDphi + theDeta*theDeta);
}

Double_t cleverRange(TH1* h,Float_t fac=1.2, Float_t minY=1.e-3)
{
   Float_t maxY =  fac * h->GetBinContent(h->GetMaximumBin());
   //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
   h->SetAxisRange(minY,maxY,"Y");
   return maxY;
}


Double_t getCleverRange(TH1* h)
{
  Double_t maxY = -1000000;
  for ( Int_t ibin = 1 ; ibin <= h->GetNbinsX() ; ibin++) {
    if (maxY < h->GetBinContent(ibin) )
      maxY = h->GetBinContent(ibin);
  }
  return maxY;
}

Double_t cleverRange(TH1* h,TH1* h2, Float_t fac=1.2, Float_t minY=1.e-3)
{
  Float_t maxY1 =  fac * h->GetBinContent(h->GetMaximumBin());
  Float_t maxY2 =  fac * h2->GetBinContent(h2->GetMaximumBin());

  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
  h->SetAxisRange(minY,max(maxY1,maxY2),"Y");
  h2->SetAxisRange(minY,max(maxY1,maxY2),"Y");
  return max(maxY1,maxY2);
}

#endif
