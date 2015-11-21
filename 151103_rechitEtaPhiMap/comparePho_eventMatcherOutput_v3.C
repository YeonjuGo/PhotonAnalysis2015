// compare photon variables on between 53X and 75X CMSSW version.
// v3 : spike rejected (only for barrel)
//    : hcalIso added
//    : argument ptThr added
//    : modified to compare data and MC on 75X
//
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

void compareTwo(TH1* h1=0, TH1* h2=0,double xmin=0.0, double xmax=200.0, const char* name="_pt15_etaBarrel",double ymax=-1);
void comparePho_eventMatcherOutput_v3(int ptThr = 15, bool isGED=0)
{
    TH1::SetDefaultSumw2();
    gStyle->SetLabelSize(0.03,"Y");
    gStyle -> SetTitleYSize(0.05);
    gStyle -> SetTitleXSize(0.05);
    const char* fname[2];
    fname[0] = "/mnt/hadoop/cms/store/user/tatar/HIHighPt/HiForest_HIHighPt_photon30_HIRun2011-v1.root"; // 75X data 
    fname[1] = "/mnt/hadoop/cms/store/user/tatar/Pyquen_Unquenched_AllQCDPhoton30_PhotonFilter20GeV_eta24_TuneZ2_PbPb_5020GeV/Pyquen_Unquenched_AllQCDPhoton30_FOREST_753p1.root"; // 75X MC
    TFile *inf[2];
    TTree* tpho[2];
    const char* stGED = "";
    if(isGED) stGED = "GED";
    for(int i=0;i<2;i++){
        inf[i] = new TFile(fname[i],"READ");
        tpho[i] = (TTree*)inf[i]->Get(Form("ggHiNtuplizer%s/EventTree",stGED));
    }
    
    const int Neta = 4; 
    TCut ptCut;
    TCut etaCut[Neta];
    TCut spikeCut;
    TCut totCut[Neta];//[eta]
    const char* tCut;
    ptCut = Form("phoEt>%d",ptThr);
    spikeCut = "pho_swissCrx<0.9 && abs(pho_seedTime)<3";
    etaCut[0] = "abs(phoEta)<1.44";
    etaCut[1] = "abs(phoEta)>1.44 && abs(phoEta)<2.0";
    etaCut[2] = "abs(phoEta)>2.0 && abs(phoEta)<2.5";
    etaCut[3] = "abs(phoEta)>=0";
    for(int ieta=0; ieta<Neta;ieta++){
        // adjust the spike cut on only the barrel. 
        if(ieta==0) totCut[ieta] = ptCut && etaCut[ieta] && spikeCut;
        else totCut[ieta] = ptCut && etaCut[ieta];
    }   

    const double xMin(0.0), xMax(2000.0);
    const int xBin=50;
    double energyBin[] = {0,20,40,60,80,100,120,140,160,180,200,220,240,300,500,700,900,1200,1600,2000};
    int Nenergy = sizeof(energyBin)/sizeof(double)-1;
    double ptBin[] = {0,5,10,15,20,25,30,35,40,45,50,60,70,85,110,120,140,160,200};
    int Npt = sizeof(ptBin)/sizeof(double)-1;

    TH1D* phopt[2][Neta];
    TH1D* phoE[2][Neta];
    TH1D* phoRawE[2][Neta];
    TH1D* phoEoverRawE[2][Neta];
    TH1D* phoSigmaIetaIeta[2][Neta];
    TH1D* phoEtaWidth[2][Neta];
    TH1D* phoPhiWidth[2][Neta];
    TH1D* phoBrem[2][Neta];
    TH1D* phoR9[2][Neta];
    TH1D* phoHoverE[2][Neta];
    TH1D* phoHcalIso[2][Neta];
    cout << "ss"<< endl;

    for(int i=0; i<2;i++){
        for(int ieta=0; ieta<Neta;ieta++){
            phopt[i][ieta] = new TH1D(Form("phopt_%d_ieta%d",i,ieta),";Photon p_{T} (GeV);",Npt,ptBin);
            phoE[i][ieta]= new TH1D(Form("phoE_%d_ieta%d",i,ieta),";Photon Energy (GeV);",Nenergy,energyBin);
            phoRawE[i][ieta] = new TH1D(Form("phoRawE_%d_ieta%d",i,ieta),";Photon Raw Energy (GeV);",Nenergy,energyBin);
            phoEoverRawE[i][ieta] = new TH1D(Form("phoEOverRawE_%d_ieta%d",i,ieta),";Photon Energy/RawEnergy;",xBin,0.95,1.2);
            phoSigmaIetaIeta[i][ieta] = new TH1D(Form("phoSigmaIetaIeta_%d_ieta%d",i,ieta),";Photon #sigma_{I#etaI#eta};",xBin,0.0,0.1);
            phoEtaWidth[i][ieta] = new TH1D(Form("phoEtaWidth_%d_ieta%d",i,ieta),";Photon SC Eta Width;",xBin,0.0,0.27);
            phoPhiWidth[i][ieta] = new TH1D(Form("phoPhiWidth_%d_ieta%d",i,ieta),";Photon SC Phi Width;",xBin,0.0,0.35);
            phoBrem[i][ieta] = new TH1D(Form("phoBrem_%d_ieta%d",i,ieta),";Photon SC PhiWidth/EtaWidth;",xBin,0.0,0.1);
            phoR9[i][ieta] = new TH1D(Form("phoR9_%d_ieta%d",i,ieta),";Photon R9;",xBin,0.0,1.0);
            phoHoverE[i][ieta] = new TH1D(Form("phoHoverE_%d_ieta%d",i,ieta),";Photon H/E;",xBin,0.0,1.0);
            phoHcalIso[i][ieta] = new TH1D(Form("phoHcalIso_%d_ieta%d",i,ieta),";Photon hcalIso;",xBin,-30,100);
        }
    }

    for(int i=0; i<2;i++){
        for(int ieta=0; ieta<Neta;ieta++){
            ///////////////////////////////////////////////////
            ///////////////////////////////////////////////////
            // Fill the histograms
            tpho[i]->Draw(Form("phoEt>>+%s",phopt[i][ieta]->GetName()),totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoE>>+%s",phoE[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoSCRawE>>%s",phoRawE[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoE/phoSCRawE>>%s",phoEoverRawE[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoSigmaIEtaIEta>>%s",phoSigmaIetaIeta[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoSCEtaWidth>>%s",phoEtaWidth[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoSCPhiWidth>>%s",phoPhiWidth[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoSCPhiWidth/phoSCEtaWidth>>%s",phoBrem[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoR9>>%s",phoR9[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("phoHoverE>>%s",phoHoverE[i][ieta]->GetName()), totCut[ieta].GetTitle());
            tpho[i]->Draw(Form("pho_hcalRechitIsoR4>>%s",phoHcalIso[i][ieta]->GetName()), totCut[ieta].GetTitle());
        }
    }
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    // DRAW using compareTwo function 

    for(int ieta=0; ieta<Neta;ieta++){
        const char* Name = Form("_pt%d_spikeRejected_Barrel",ptThr);
        if(ieta==1) Name = Form("_pt%d_spikeRejected_endcapEtaTo2",ptThr);
        if(ieta==2) Name = Form("_pt%d_spikeRejected_endcapEtaOver2To2p5",ptThr);
        if(ieta==3) Name = Form("_pt%d_spikeRejected_totalEta",ptThr);


        compareTwo(phopt[0][ieta], phopt[1][ieta],0,200,Name);
        compareTwo(phoE[0][ieta], phoE[1][ieta],xMin,energyBin[Nenergy],Name);
        compareTwo(phoRawE[0][ieta], phoRawE[1][ieta],xMin,energyBin[Nenergy],Name);
        compareTwo(phoEoverRawE[0][ieta], phoEoverRawE[1][ieta],0.95,1.2,Name);
        compareTwo(phoSigmaIetaIeta[0][ieta], phoSigmaIetaIeta[1][ieta],0.0,0.1,Name);
        compareTwo(phoEtaWidth[0][ieta], phoEtaWidth[1][ieta],0.0,0.27,Name);
        compareTwo(phoPhiWidth[0][ieta], phoPhiWidth[1][ieta],0.0,0.35,Name);
        compareTwo(phoR9[0][ieta], phoR9[1][ieta],0,1.0,Name);
        compareTwo(phoHoverE[0][ieta], phoHoverE[1][ieta],0.0,1.8,Name);
        compareTwo(phoHcalIso[0][ieta], phoHcalIso[1][ieta],-30.0,100,Name);
    }

}

void compareTwo(TH1* h1, TH1* h2, double xmin, double xmax, const char* name, double ymax)  {
    TCanvas* c=  new TCanvas(Form("c_%s",h1->GetName()),"", 400,800);
    c->Divide(1,2);
    c->cd(1);
    h1->Sumw2();
    h2->Sumw2();
    h1->Scale( 1. / h1->Integral());
    h2->Scale( 1. / h2->Integral());

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
    h1->Divide(h2);
    if(ymax!=-1) h1->SetAxisRange(0,ymax,"Y");
    else {
        h1->SetMinimum(0.0);
        if(h1->GetMaximum()>20) h1->SetMaximum(20);
    }
    h1->SetYTitle("DATA/MC Ratio");
    h1->DrawCopy();
    jumSun(0,1,xmax,1);
    c -> SaveAs(Form("pdf_DATAandMC/%s%s.pdf",h1->GetName(),name)); 
}

