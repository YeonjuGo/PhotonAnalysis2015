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

void prettyPlots()
{
  gStyle->SetOptStat(0);

  const int NTRIGCalo = 5;
  const int NTRIGPF = 4;
  TGraphAsymmErrors *ak4CaloTurnOns[NTRIGCalo];
  TGraphAsymmErrors *ak4PFTurnOns[NTRIGPF];

  TString ak4CaloName[NTRIGCalo] = {"HLT_AK4CaloJet30_v2", "HLT_AK4CaloJet40_v2", "HLT_AK4CaloJet50_v2", "HLT_AK4CaloJet80_v2", "HLT_AK4CaloJet100_v2"};
  TString ak4PFName[NTRIGPF] = {"HLT_AK4PFJet30_v2", "HLT_AK4PFJet50_v2", "HLT_AK4PFJet80_v2", "HLT_AK4PFJet100_v2"};

  Int_t trigColors[NTRIGCalo] = {1, kBlue, kRed, 90, kMagenta};

  TFile *inFile = TFile::Open("jetTurnOn.root");

  for(int i = 0; i < NTRIGCalo; i++){
    ak4CaloTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak4CaloName[i]+"_asymm");
    ak4CaloTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak4CaloTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < NTRIGPF; i++){
    ak4PFTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak4PFName[i]+"_asymm");
    ak4PFTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak4PFTurnOns[i]->SetLineColor(trigColors[i]);
  }

  TH1D *hEmpty = new TH1D("hEmpty",";p_{T}^{reco};Efficiency",20,0,150);
  hEmpty->SetMaximum(1.2);
  hEmpty->SetMinimum(0.0);

  TCanvas *c1[2];
  TString canvName[2] = {"ak4CaloTrig", "ak4PFTrig"};
  TLine *oneLine = new TLine(0,1,150,1);
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
  for(int i = 0; i < NTRIGCalo; ++i){
    //    ak4CaloTurnOns[i]->UseCurrentStyle(0);
    ak4CaloTurnOns[i]->Draw("p e");
    leg[0]->AddEntry(ak4CaloTurnOns[i], ak4CaloName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 5.0");


  leg[0]->Draw("SAME");

  c1[1]->cd();
  for(int i = 0; i < NTRIGPF; ++i){
    //    ak4PFTurnOns[i]->UseCurrentStyle(1);
    ak4PFTurnOns[i]->Draw("p e");
    leg[1]->AddEntry(ak4PFTurnOns[i], ak4PFName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 5.0");
  leg[1]->Draw("SAME");

  TFile* outFile_p = new TFile("jetTurnOnPlots.root", "UPDATE");
  for(Int_t iter = 0; iter < 2; iter++){
    c1[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(c1[iter], Form("pdfDir/RECOETA5_NEUTRINO%s", canvName[iter].Data()), "pdf");
  }
  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < 2; iter++){
    delete leg[iter];
  }

  delete oneLine;

  for(Int_t iter = 0; iter < 2; iter++){
    delete c1[iter];
  }
  delete hEmpty;


  inFile->Close();

  return;
}
