#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLatex.h"
#include "stdio.h"
#include "../yjUtility.h"
#include "../HIUtils/histoUtil.h"
//last forward run is 211256


//const Double_t hfBins[] = {0, 20, 30, 1000}; //last entry is upper bound on last bin
//const Int_t nhfBins = 3;

//int returnHFBin(double hf);
void compareTwo(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, const char* theCut="");
void compareRechit()
{
    TH1::SetDefaultSumw2();
    const char* fname[2];
    fname[0] = "/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_AllQCDPhoton30_ecalGlobal_755p1.root";
    fname[1] = "/afs/cern.ch/work/y/ygo/public/ecalLocalReco/forest_AllQCDPhoton30_ecalMultifit_755p1.root";
    TFile* inf[2];
    TTree* tee[2];
    TTree* teb[2];
    //TTree* tevt[2];
    //TTree* thbhe[2];
    TTree* tpho[2];
    const char* ver = "GlobalReco";
    for(int i=0;i<2;i++){
        if(i==1) ver = "MultiFit";
        inf[i] = new TFile(fname[i], "READ");
        tee[i] = (TTree*)inf[i]->Get(Form("rechitanalyzer/ee"));
        teb[i] = (TTree*)inf[i]->Get(Form("rechitanalyzer/eb"));
        //thbhe[i] = (TTree*)inf[i]->Get(Form("rechitAnalyzer/hbhe"));
        //tevt[i] = (TTree*)inf[i]->Get("hiEvtAnalyzer/HiTree");
        tpho[i] = (TTree*)inf[i]->Get(Form("ggHiNtuplizer/EventTree"));
    }
    compareTwo(tee[0],tee[1],"e",50,0,400);
    compareTwo(tee[0],tee[1],"et",50,0,100);
    compareTwo(tee[0],tee[1],"phi",50,-3.14,3.14);
    compareTwo(tee[0],tee[1],"eta",50,-3,3);
    compareTwo(teb[0],teb[1],"e",50,0,400,"");
    compareTwo(teb[0],teb[1],"et",50,0,100);
    compareTwo(teb[0],teb[1],"phi",50,-3.14,3.14);
    compareTwo(teb[0],teb[1],"eta",50,-3,3);

}

void compareTwo(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, const char* theCut)  {
    TCanvas* cc=  new TCanvas(Form("c_%s_%s",var.Data(),t1->GetName()),"", 400,800);
    cc->Divide(1,2);
    cc->cd(1);
    TH1D* h1 = new TH1D(Form("h1_%s_%s",var.Data(),t1->GetName()), Form(";%s;",var.Data()), nBins,xMin,xMax);
    TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%s",var.Data(),t1->GetName()));
    h1->Sumw2();
    h2->Sumw2();
    t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), theCut);
    t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), theCut);
    //  h1->Scale( 1. / h1->Integral());
    //  h2->Scale( 1. / h2->Integral());

    // default plotting options
    h1->SetMarkerStyle(kFullSquare);
    h2->SetMarkerStyle(kFullCircle);

    h1->SetMarkerColor(kBlack);
    h2->SetMarkerColor(kRed);

    h2->SetMarkerSize(h1->GetMarkerSize()*0.8);     // to distinguish between points when they overlap

    h1->SetLineColor(4);
    h1->SetMarkerColor(4);
    h2->SetLineColor(2);
    h2->SetMarkerColor(2);

    // set Y-axis ranges
    // default Y-axis range is that of h1
    // make sure that both plots will not run out of y-axis
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    double min1 = h1->GetMinimum();
    double min2 = h2->GetMinimum();

    if (max2 > max1)
        h1->SetMaximum(max2+TMath::Abs(max2)*0.2);

    if (min2 < min1)
        h1->SetMinimum(min2-TMath::Abs(min2)*0.2);

    h2->DrawCopy();
    h1->DrawCopy("hist same");
    TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
    //legStyle(leg);
    leg->AddEntry(h1,Form("Global"), "l");
    leg->AddEntry(h2,Form("Multifit"));
    leg->Draw();
    cc->cd(2);
    h2->SetAxisRange(0.5,1.5,"Y");
    h2->Divide(h1);
    h2->SetYTitle("GlobalReco/MultiFit Ratio");
    h2->DrawCopy();
    jumSun(xMin,1,xMax,1);
    cc-> SaveAs(Form("pdf/%s.pdf",cc->GetName())); 
    //return c; 
}

int main(){
    compareRechit();
    return 0;
}
