#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "TLatex.h"
#include "stdio.h"
#include "yjUtility.h"

void compareTwo(TFile* f1, TFile* f2, const char* tree="rechitanalyzer/tower", TString var="et", int nBins=40, double xMin=0, double xMax=20, TCut theCut1="", TCut theCut2="", const char* cap="")  {
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0); 
    TTree* t1 = (TTree*) f1->Get(tree);
    TTree* t2 = (TTree*) f2->Get(tree);
    TCanvas* c=  new TCanvas(Form("c_%s%s",var.Data(),cap),"", 400,800);
	c->Divide(1,2);
	c->cd(1);
	gPad->SetLogy();
	TH1D* h1 = new TH1D(Form("h1_%s%s",var.Data(),cap), Form(";%s;",var.Data()), nBins,xMin,xMax);
	TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s%s",var.Data(),cap));
	h1->Sumw2();
	h2->Sumw2();
    cout << "bb" << endl;
	t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut1);
    cout << "bb" << endl;
	t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut2);	
    cout << "bb" << endl;
	h1->Scale( 1. / t1->GetEntries(theCut1));
	h2->Scale( 1. / t2->GetEntries(theCut2));
    cout << "bb" << endl;
    SetHistColor(h1,2);
    SetHistColor(h2,1);
    cout << "bb" << endl;
    h1->SetMarkerStyle(20);
	h1->SetMarkerSize(0.8);
    cout << "bb" << endl;
	double range = cleverRange(h1,h2,1.5,1.e-4);
	h1->DrawCopy("L");
    cout << "bb" << endl;
	h2->DrawCopy("hist same");
	c->cd(2);
    cout << "bb" << endl;
	h2->Divide(h1);
	h2->SetYTitle("Ratio");
    cout << "bb" << endl;
	h2->GetYaxis()->SetRangeUser(0.0,2.0);
    SetHistColor(h2,2);
    cout << "bb" << endl;
    h2->SetMarkerStyle(20);
	h2->SetMarkerSize(0.8);
	h2->DrawCopy("");
	//h2->DrawCopy("le1");
	jumSun(xMin,1,xMax,1);
	c-> SaveAs(Form("pdf/compare_%s%s.pdf",var.Data(),cap));
}
