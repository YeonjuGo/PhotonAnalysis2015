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
#include "TCut.h"
#include "stdio.h"
#include "../yjUtility.h"
//#include "../HIUtils/histoUtil.h"
//last forward run is 211256


//const Double_t hfBins[] = {0, 20, 30, 1000}; //last entry is upper bound on last bin
//const Int_t nhfBins = 3;

//int returnHFBin(double hf);

void compareTwo(TH1* h1=0, TH1* h2=0,double xmin=0.0, double xmax=200.0, const char* name="_pt15_etaBarrel",double ymax=3.0);
void comparePho_eventMatcherOutput()
{
    TH1::SetDefaultSumw2();
    gStyle->SetLabelSize(0.03,"Y");
    gStyle -> SetTitleYSize(0.05);
    gStyle -> SetTitleXSize(0.06);
    TFile *inf1 = new TFile("/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/151103_rechitEtaPhiMap/skimFiles/eventMatched_53X_75X_biransSample.root");
    //TFile *inf1 = new TFile("/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/151103_rechitEtaPhiMap/skimFiles/eventMatched_53X_75X.root");
    TTree* tevt = (TTree*)inf1->Get("HiTree");
    TTree* tpho[2];
    const char* ver = "53";
    for(int i=0;i<2;i++){
        if(i==1) ver = "75";
        tpho[i] = (TTree*)inf1->Get(Form("pho%sx",ver));
    }

    const double xMin(0.0), xMax(2000.0);
    const int xBin=50;

    double energyBin[] = {0,20,40,60,80,100,120,140,160,180,200,220,240,300,500,700,900,1200,1600,2000};
    int Nenergy = sizeof(energyBin)/sizeof(double)-1;
    double ptBin[] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,140,160,180,200};
    int Npt = sizeof(ptBin)/sizeof(double)-1;

    const int Neta = 4;
    TCut ptCut[2];
    TCut etaCut[2][Neta];
    TCut spikeCut[2];
    TCut totCut[2][Neta];//[sample][eta]
    ptCut[0] = "pt>15";
    ptCut[1] = "phoEt>15";
    spikeCut[0] = "swissCrx<0.9 && abs(seedTime)<3";
    spikeCut[1] = "pho_swissCrx<0.9 && abs(pho_seedTime)<3";
    etaCut[0][0] = "pt>15 && abs(eta)<1.44";
    etaCut[1][0] = "phoEt>15 && abs(phoEta)<1.44";
    etaCut[0][1] = "pt>15 && abs(eta)>1.44 && abs(eta)<2.0";
    etaCut[1][1] = "phoEt>15 && abs(phoEta)>1.44 && abs(phoEta<2.0)";
    etaCut[0][2] = "pt>15 && abs(eta)>2.0 && abs(eta)<2.5";
    etaCut[1][2] = "phoEt>15 && abs(phoEta)>2.0 && abs(phoEta)<2.5";
    etaCut[0][3] = "pt>15 && abs(eta)>=0";
    etaCut[1][3] = "phoEt>15 && abs(phoEta)>=0";

    for(int i=0; i<2;i++){
        for(int ieta=0; ieta<Neta;ieta++){
            totCut[i][ieta] = etaCut[i][ieta] && spikeCut[i];
        }
    }

    TH1D* phopt[2][Neta];
    TH1D* phoE[2][Neta];
    TH1D* phoRawE[2][Neta];
    TH1D* phoCoshRawE[2][Neta];
    TH1D* phoEoverRawE[2][Neta];
    TH1D* phoSigmaIetaIeta[2][Neta];
    TH1D* phoEtaWidth[2][Neta];
    TH1D* phoPhiWidth[2][Neta];
    TH1D* phoBrem[2][Neta];
    TH1D* phoR9[2][Neta];
    TH1D* phoHoverE[2][Neta];
    cout << "ss"<< endl;

    for(int i=0; i<2;i++){
        for(int ieta=0; ieta<Neta;ieta++){
            phopt[i][ieta] = new TH1D(Form("phopt_%d_ieta%d",i,ieta),";Photon p_{T} (GeV);",Npt,ptBin);
            phoE[i][ieta]= new TH1D(Form("phoE_%d_ieta%d",i,ieta),";Photon Energy (GeV);",Nenergy,energyBin);
            phoRawE[i][ieta] = new TH1D(Form("phoRawE_%d_ieta%d",i,ieta),";Photon Raw Energy (GeV);",Nenergy,energyBin);
            //phoCoshRawE[i][ieta] = new TH1D(Form("phoCoshRawE_%d_ieta%d",i,ieta),";Photon Raw Energy (GeV);",xBin,xMin,xMax);
            phoEoverRawE[i][ieta] = new TH1D(Form("phoEOverRawE_%d_ieta%d",i,ieta),";Photon Energy/RawEnergy;",xBin,0.95,1.2);
            phoSigmaIetaIeta[i][ieta] = new TH1D(Form("phoSigmaIetaIeta_%d_ieta%d",i,ieta),";Photon #sigma_{I#etaI#eta};",xBin,0.0,0.1);
            phoEtaWidth[i][ieta] = new TH1D(Form("phoEtaWidth_%d_ieta%d",i,ieta),";Photon SC Eta Width;",xBin,0.0,0.27);
            phoPhiWidth[i][ieta] = new TH1D(Form("phoPhiWidth_%d_ieta%d",i,ieta),";Photon SC Phi Width;",xBin,0.0,0.35);
            phoBrem[i][ieta] = new TH1D(Form("phoBrem_%d_ieta%d",i,ieta),";Photon SC PhiWidth/EtaWidth;",xBin,0.0,0.1);
            phoR9[i][ieta] = new TH1D(Form("phoR9_%d_ieta%d",i,ieta),";Photon R9;",xBin,0.0,1.0);
            phoHoverE[i][ieta] = new TH1D(Form("phoHoverE_%d_ieta%d",i,ieta),";Photon H/E;",xBin,0.0,1.0);
        }
    }
    
    cout << "aa"<< endl;
    for(int ieta=0; ieta<Neta;ieta++){
        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
        // Fill the histograms 
        tpho[0]->Draw(Form("pt>>%s",phopt[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoEt>>%s",phopt[1][ieta]->GetName()),totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("energy>>%s",phoE[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoE>>%s",phoE[1][ieta]->GetName()), totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("rawEnergy>>%s",phoRawE[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoSCRawE>>%s",phoRawE[1][ieta]->GetName()), totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("energy/rawEnergy>>%s",phoEoverRawE[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoE/phoSCRawE>>%s",phoEoverRawE[1][ieta]->GetName()), totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("sigmaIetaIeta>>%s",phoSigmaIetaIeta[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoSigmaIEtaIEta>>%s",phoSigmaIetaIeta[1][ieta]->GetName()), totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("etaWidth>>%s",phoEtaWidth[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoSCEtaWidth>>%s",phoEtaWidth[1][ieta]->GetName()), totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("phiWidth>>%s",phoPhiWidth[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoSCPhiWidth>>%s",phoPhiWidth[1][ieta]->GetName()), totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("phiWidth/etaWidth>>%s",phoBrem[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoSCPhiWidth/phoSCEtaWidth>>%s",phoBrem[1][ieta]->GetName()), totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("r9>>%s",phoR9[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoR9>>%s",phoR9[1][ieta]->GetName()), totCut[1][ieta].GetTitle());
        tpho[0]->Draw(Form("hadronicOverEm>>%s",phoHoverE[0][ieta]->GetName()), totCut[0][ieta].GetTitle());
        tpho[1]->Draw(Form("phoHoverE>>%s",phoHoverE[1][ieta]->GetName()), totCut[1][ieta].GetTitle());

        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
        // DRAW using compareTwo function 

        const char* Name = "Barrel";
        if(ieta==1) Name = "endcapEtaTo2";
        if(ieta==2) Name = "endcapEtaOver2To2p5";
        if(ieta==3) Name = "totalEta";


        compareTwo(phopt[0][ieta], phopt[1][ieta],0,200,Form("_pt15_%s",Name));
        if(ieta==1) compareTwo(phoE[0][ieta], phoE[1][ieta],xMin,xMax,Form("_pt15_%s",Name),10);
        else compareTwo(phoE[0][ieta], phoE[1][ieta],xMin,xMax,Form("_pt15_%s",Name));
        if(ieta==1) compareTwo(phoRawE[0][ieta], phoRawE[1][ieta],xMin,xMax,Form("_pt15_%s",Name),10);
        else compareTwo(phoRawE[0][ieta], phoRawE[1][ieta],xMin,xMax,Form("_pt15_%s",Name));
        //compareTwo(phoCoshRawE[0][ieta], phoCoshRawE[1][ieta],xMin,xMax,Form("_pt15_%s",Name));
        compareTwo(phoEoverRawE[0][ieta], phoEoverRawE[1][ieta],0.95,1.2,Form("_pt15_%s",Name));
        compareTwo(phoSigmaIetaIeta[0][ieta], phoSigmaIetaIeta[1][ieta],0.0,0.1,Form("_pt15_%s",Name),9);
        compareTwo(phoEtaWidth[0][ieta], phoEtaWidth[1][ieta],0.0,0.27,Form("_pt15_%s",Name));
        compareTwo(phoPhiWidth[0][ieta], phoPhiWidth[1][ieta],0.0,0.35,Form("_pt15_%s",Name));
        compareTwo(phoR9[0][ieta], phoR9[1][ieta],0,1.0,Form("_pt15_%s",Name),30);
        compareTwo(phoHoverE[0][ieta], phoHoverE[1][ieta],0.0,1.8,Form("_pt15_%s",Name),7);
    }

}

void compareTwo(TH1* h1, TH1* h2, double xmin, double xmax, const char* name, double ymax)  {
    TCanvas* c=  new TCanvas(Form("c_%s",h1->GetName()),"", 400,800);
    c->Divide(1,2);
    c->cd(1);
    h1->Sumw2();
    h2->Sumw2();
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

    h1->DrawCopy("hist");
    h2->DrawCopy("hist same");
    c->cd(2);
    h2->SetAxisRange(0,ymax,"Y");
    h2->Divide(h1);
    h2->SetYTitle("75X/53X Ratio");
    h2->DrawCopy();
    jumSun(xmin,1,xmax,1);
    c -> SaveAs(Form("pdf/%s%s.pdf",h1->GetName(),name)); 
}

