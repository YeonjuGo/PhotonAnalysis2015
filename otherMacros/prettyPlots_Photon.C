#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TDatime.h>
#include <TLatex.h>

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif") {
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
}

void prettyPlots_Photon()
{
  gStyle->SetOptStat(0);

  const int nL1Trig = 9;
  TGraphAsymmErrors *L1TurnOns[nL1Trig];

  const int nHLTTrig = 5;
  TGraphAsymmErrors *photonTurnOns[nHLTTrig];

  TString L1Name[nL1Trig] = {"L1_SingleEG2", "L1_SingleEG5", "L1_SingleEG10", "L1_SingleEG15", "L1_SingleEG20", "L1_SingleEG25", "L1_SingleEG30", "L1_SingleEG35", "L1_SingleEG40"};

  TString photonName[nHLTTrig] = {"HLT_HISinglePhoton10_v2", "HLT_HISinglePhoton15_v2", "HLT_HISinglePhoton20_v2", "HLT_HISinglePhoton40_v2", "HLT_HISinglePhoton60_v2"};


  Int_t trigHLTColors[nHLTTrig] = {1, kBlue, kRed, 90, kMagenta};
  Int_t trigL1Colors[nL1Trig] = {1, 1, kBlue, kBlue, kRed, kRed, 90, 90, kMagenta};

  TFile *inFile = TFile::Open("photonTurnOn.root");

  for(int i = 0; i < nL1Trig; i++){
    L1TurnOns[i] = (TGraphAsymmErrors *)inFile->Get(L1Name[i]+"_asymm");
    L1TurnOns[i]->SetMarkerColor(trigL1Colors[i]);
    L1TurnOns[i]->SetLineColor(trigL1Colors[i]);
  }

  for(int i = 0; i < nHLTTrig; i++){
    photonTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(photonName[i]+"_asymm");
    photonTurnOns[i]->SetMarkerColor(trigHLTColors[i]);
    photonTurnOns[i]->SetLineColor(trigHLTColors[i]);
  }

  TH1D *hEmpty = new TH1D("hEmpty",";p_{T}^{reco};Efficiency",20,0,80);
  hEmpty->SetMaximum(1.2);
  hEmpty->SetMinimum(0.0);

  TCanvas *c1[2];
  TString canvName[2] = {"L1Trig", "photonTrig"};
  TLine *oneLine = new TLine(0,1,80,1);
  oneLine->SetLineStyle(2);
  TLegend *leg[2];

  for(Int_t iter = 0; iter < 2; iter++){
    c1[iter] = new TCanvas(canvName[iter], canvName[iter], 700, 700);//("c1","c1",700,700);
    
    hEmpty->DrawCopy();
    oneLine->Draw();

    leg[iter] = new TLegend(0.65, 0.20, 0.98, 0.45);
    leg[iter]->SetFillColor(0);
    leg[iter]->SetTextFont(43);
    leg[iter]->SetTextSize(16);
  }
  
  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(23);

  c1[0]->cd();  

  for(int i = 0; i < nL1Trig; ++i){
    if(i%2 != 0) continue;
    L1TurnOns[i]->Draw("p e");
    leg[0]->AddEntry(L1TurnOns[i], L1Name[i], "p l");
  }

  label_p->DrawLatex(.75, .50, "|#eta| < 1.44");
  leg[0]->Draw("SAME");

  c1[1]->cd();

  for(int i = 0; i < nHLTTrig; ++i){
    //    photonTurnOns[i]->UseCurrentStyle(0);
    photonTurnOns[i]->Draw("p e");
    leg[1]->AddEntry(photonTurnOns[i], photonName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 1.44");

  leg[1]->Draw("SAME");

  TFile* outFile_p = new TFile("photonTurnOnPlots.root", "UPDATE");
  c1[0]->Write("", TObject::kOverwrite);
  claverCanvasSaving(c1[0], Form("pdfDir/RECOETA144_NEUTRINO%s", canvName[0].Data()), "pdf");

  c1[1]->Write("", TObject::kOverwrite);
  claverCanvasSaving(c1[1], Form("pdfDir/RECOETA144_NEUTRINO%s", canvName[1].Data()), "pdf");
  outFile_p->Close();
  delete outFile_p;

  delete leg[0];
  delete leg[1];

  delete oneLine;

  delete c1[0];
  delete c1[1];
  
  delete hEmpty;


  inFile->Close();

  return;
}
