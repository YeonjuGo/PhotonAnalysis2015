#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TDatime.h>
#include <TLatex.h>

#include <iostream>

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif"){
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
  return;
}

void prettyPlots()
{
  gStyle->SetOptStat(0);

  const int n2Trig = 5;
  const int n3Trig = 5;
  const int n4Trig = 5;
  const int n5Trig = 5;
  const int nPhotonTrig = 5;

  TH1F* rateHists_p[4];

  TGraphAsymmErrors *ak2CaloTurnOns[n2Trig];
  TGraphAsymmErrors *ak3CaloTurnOns[n3Trig];
  TGraphAsymmErrors *ak4CaloTurnOns[n4Trig];
  TGraphAsymmErrors *ak5CaloTurnOns[n5Trig];

  TGraphAsymmErrors *ak2CaloTrkTurnOns[n2Trig];
  TGraphAsymmErrors *ak3CaloTrkTurnOns[n3Trig];
  TGraphAsymmErrors *ak4CaloTrkTurnOns[n4Trig];
  TGraphAsymmErrors *ak5CaloTrkTurnOns[n5Trig];

  TGraphAsymmErrors *ak2CaloGenTurnOns[n2Trig];
  TGraphAsymmErrors *ak3CaloGenTurnOns[n3Trig];
  TGraphAsymmErrors *ak4CaloGenTurnOns[n4Trig];
  TGraphAsymmErrors *ak5CaloGenTurnOns[n5Trig];

  TGraphAsymmErrors *photonTurnOns[nPhotonTrig];

  TString ak2CaloName[n2Trig] = {"HLT_PuAK2CaloJet40_v1", "HLT_PuAK2CaloJet60_v1", "HLT_PuAK2CaloJet80_v1", "HLT_PuAK2CaloJet100_v1", "HLT_PuAK2CaloJet120_v1"};
  TString ak3CaloName[n3Trig] = {"HLT_PuAK3CaloJet40_v1", "HLT_PuAK3CaloJet60_v1", "HLT_PuAK3CaloJet80_v1", "HLT_PuAK3CaloJet100_v1", "HLT_PuAK3CaloJet120_v1"};
  TString ak4CaloName[n4Trig] = {"HLT_PuAK4CaloJet40_v1", "HLT_PuAK4CaloJet60_v1", "HLT_PuAK4CaloJet80_v1", "HLT_PuAK4CaloJet100_v1", "HLT_PuAK4CaloJet120_v1"};
  TString ak5CaloName[n5Trig] = {"HLT_PuAK5CaloJet40_v1", "HLT_PuAK5CaloJet60_v1", "HLT_PuAK5CaloJet80_v1", "HLT_PuAK5CaloJet100_v1", "HLT_PuAK5CaloJet120_v1"};

  TString ak2CaloTrkName[n2Trig] = {"HLT_PuAK2CaloTrkJet40_v1", "HLT_PuAK2CaloTrkJet60_v1", "HLT_PuAK2CaloTrkJet80_v1", "HLT_PuAK2CaloTrkJet100_v1", "HLT_PuAK2CaloTrkJet120_v1"};
  TString ak3CaloTrkName[n3Trig] = {"HLT_PuAK3CaloTrkJet40_v1", "HLT_PuAK3CaloTrkJet60_v1", "HLT_PuAK3CaloTrkJet80_v1", "HLT_PuAK3CaloTrkJet100_v1", "HLT_PuAK3CaloTrkJet120_v1"};
  TString ak4CaloTrkName[n4Trig] = {"HLT_PuAK4CaloTrkJet40_v1", "HLT_PuAK4CaloTrkJet60_v1", "HLT_PuAK4CaloTrkJet80_v1", "HLT_PuAK4CaloTrkJet100_v1", "HLT_PuAK4CaloTrkJet120_v1"};
  TString ak5CaloTrkName[n5Trig] = {"HLT_PuAK5CaloTrkJet40_v1", "HLT_PuAK5CaloTrkJet60_v1", "HLT_PuAK5CaloTrkJet80_v1", "HLT_PuAK5CaloTrkJet100_v1", "HLT_PuAK5CaloTrkJet120_v1"};

  TString ak2CaloGenName[n2Trig] = {"HLT_PuAK2CaloGenJet40_v1", "HLT_PuAK2CaloGenJet60_v1", "HLT_PuAK2CaloGenJet80_v1", "HLT_PuAK2CaloGenJet100_v1", "HLT_PuAK2CaloGenJet120_v1"};
  TString ak3CaloGenName[n3Trig] = {"HLT_PuAK3CaloGenJet40_v1", "HLT_PuAK3CaloGenJet60_v1", "HLT_PuAK3CaloGenJet80_v1", "HLT_PuAK3CaloGenJet100_v1", "HLT_PuAK3CaloGenJet120_v1"};
  TString ak4CaloGenName[n4Trig] = {"HLT_PuAK4CaloGenJet40_v1", "HLT_PuAK4CaloGenJet60_v1", "HLT_PuAK4CaloGenJet80_v1", "HLT_PuAK4CaloGenJet100_v1", "HLT_PuAK4CaloGenJet120_v1"};
  TString ak5CaloGenName[n5Trig] = {"HLT_PuAK5CaloGenJet40_v1", "HLT_PuAK5CaloGenJet60_v1", "HLT_PuAK5CaloGenJet80_v1", "HLT_PuAK5CaloGenJet100_v1", "HLT_PuAK5CaloGenJet120_v1"};

  TString photonName[nPhotonTrig] = {"HLT_HISinglePhoton10_v1", "HLT_HISinglePhoton15_v1", "HLT_HISinglePhoton20_v1", "HLT_HISinglePhoton40_v1", "HLT_HISinglePhoton60_v1"};

  Int_t trigColors[n3Trig] = {1, kBlue, kRed, kYellow+1, kMagenta};

  TFile *inFile = TFile::Open("jetTurnOn_HI.root");

  for(Int_t iter = 0; iter < 4; iter++){
    rateHists_p[iter] = (TH1F*)inFile->Get(Form("rateHist_%d_h", iter+2));
    rateHists_p[iter]->SetMarkerColor(trigColors[iter]);
    rateHists_p[iter]->SetLineColor(trigColors[iter]);
    rateHists_p[iter]->SetMinimum(.1);
    rateHists_p[iter]->SetXTitle("Trigger p_{T}");
    rateHists_p[iter]->SetYTitle("Rate (Hz)");
    rateHists_p[iter]->GetXaxis()->CenterTitle();
    rateHists_p[iter]->GetYaxis()->CenterTitle();
  }

  for(int i = 0; i < n2Trig; i++){
    ak2CaloTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak2CaloName[i]+"_asymm");
    ak2CaloTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak2CaloTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n3Trig; i++){
    ak3CaloTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak3CaloName[i]+"_asymm");
    ak3CaloTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak3CaloTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n4Trig; i++){
    ak4CaloTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak4CaloName[i]+"_asymm");
    ak4CaloTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak4CaloTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n5Trig; i++){
    ak5CaloTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak5CaloName[i]+"_asymm");
    ak5CaloTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak5CaloTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n2Trig; i++){
    ak2CaloTrkTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak2CaloTrkName[i]+"_asymm");
    ak2CaloTrkTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak2CaloTrkTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n3Trig; i++){
    ak3CaloTrkTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak3CaloTrkName[i]+"_asymm");
    ak3CaloTrkTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak3CaloTrkTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n4Trig; i++){
    ak4CaloTrkTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak4CaloTrkName[i]+"_asymm");
    ak4CaloTrkTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak4CaloTrkTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n5Trig; i++){
    ak5CaloTrkTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak5CaloTrkName[i]+"_asymm");
    ak5CaloTrkTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak5CaloTrkTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n2Trig; i++){
    ak2CaloGenTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak2CaloGenName[i]+"_asymm");
    ak2CaloGenTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak2CaloGenTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n3Trig; i++){
    ak3CaloGenTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak3CaloGenName[i]+"_asymm");
    ak3CaloGenTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak3CaloGenTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n4Trig; i++){
    ak4CaloGenTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak4CaloGenName[i]+"_asymm");
    ak4CaloGenTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak4CaloGenTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < n5Trig; i++){
    ak5CaloGenTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(ak5CaloGenName[i]+"_asymm");
    ak5CaloGenTurnOns[i]->SetMarkerColor(trigColors[i]);
    ak5CaloGenTurnOns[i]->SetLineColor(trigColors[i]);
  }

  for(int i = 0; i < nPhotonTrig; i++){
    photonTurnOns[i] = (TGraphAsymmErrors *)inFile->Get(photonName[i]+"_asymm");
    photonTurnOns[i]->SetMarkerColor(trigColors[i]);
    photonTurnOns[i]->SetLineColor(trigColors[i]);
  }

  TH1D *hEmpty = new TH1D("hEmpty",";p_{T}^{reco};Efficiency",20, 20,160);
  hEmpty->SetMaximum(1.1);
  hEmpty->SetMinimum(0.0);
  hEmpty->SetXTitle("Jet p_{T}");
  hEmpty->GetXaxis()->SetTitleOffset(1.0);
  hEmpty->GetXaxis()->SetTitleColor(1);

  TH1D *hEmpty_T = new TH1D("hEmpty_T",";p_{T}^{reco};Efficiency",20, 0, 120);
  hEmpty_T->SetMaximum(1.1);
  hEmpty_T->SetMinimum(0.0);
  hEmpty_T->SetXTitle("Trk p_{T}");
  hEmpty_T->GetXaxis()->SetTitleOffset(1.0);
  hEmpty_T->GetXaxis()->SetTitleColor(1);

  TH1D *hEmpty_Gen = new TH1D("hEmpty_Gen",";p_{T}^{reco};Efficiency",20, 0, 120);
  hEmpty_Gen->SetMaximum(1.1);
  hEmpty_Gen->SetMinimum(0.0);
  hEmpty_Gen->SetXTitle("Gen. particle p_{T} (chg!=0)");
  hEmpty_Gen->GetXaxis()->SetTitleOffset(1.0);
  hEmpty_Gen->GetXaxis()->SetTitleColor(1);

  TH1D *hEmpty_Gam = new TH1D("hEmpty_Gam",";p_{T}^{reco};Efficiency",20,0,100);
  hEmpty_Gam->SetMaximum(1.1);
  hEmpty_Gam->SetMinimum(0.0);

  const Int_t nCanv = 13;

  TCanvas *rateCanv_p = new TCanvas("jetRateCanv_c", "jetRateCanv_c", 700, 700);
  TLegend *rateLeg_p = new TLegend(0.65, 0.60, 0.85, 0.85);
  rateLeg_p->SetFillColor(0);
  rateLeg_p->SetTextFont(43);
  rateLeg_p->SetTextSize(16);

  TLine* tenLine_p = new TLine(30, 10, 130, 10);
  TLine* oneLine_p = new TLine(30, 1, 130, 1);
  tenLine_p->SetLineStyle(2);

  const std::string drawOpt[4] = {"E1 SAME", "E1 SAME", "E1 SAME", "E1"};
  for(Int_t iter = 3; iter > -1; iter--){
    rateCanv_p->cd();
    if(iter == 3) gPad->SetLogy();

    rateHists_p[iter]->DrawCopy(drawOpt[iter].c_str());
    rateLeg_p->AddEntry(rateHists_p[iter], Form("Rate R = 0.%d", iter+2), "P L");

    if(iter == 0){
      tenLine_p->Draw();
      oneLine_p->Draw();
    }
  }

  rateLeg_p->Draw("SAME");

  TCanvas *c1[nCanv];
  TString canvName[nCanv] = {"akPu2CaloTrig", "akPu3CaloTrig", "akPu4CaloTrig", "akPu5CaloTrig", "akPu2CaloTrkTrig", "akPu3CaloTrkTrig", "akPu4CaloTrkTrig", "akPu5CaloTrkTrig", "akPu2CaloGenTrig", "akPu3CaloGenTrig", "akPu4CaloGenTrig", "akPu5CaloGenTrig", "photonTrig"};
  TLine *oneLine = new TLine(20,1,160,1);
  oneLine->SetLineStyle(2);
  TLine *oneLine_T = new TLine(0,1,120,1);
  oneLine_T->SetLineStyle(2);
  TLine *oneLine_Gen = new TLine(0,1,120,1);
  oneLine_Gen->SetLineStyle(2);
  TLine *oneLine_Gam = new TLine(0,1,100,1);
  oneLine_Gam->SetLineStyle(2);
  TLegend *leg[nCanv];

  for(Int_t iter = 0; iter < nCanv; iter++){
    c1[iter] = new TCanvas(canvName[iter], canvName[iter], 700, 700);//("c1","c1",700,700);

    if(iter < 4){
      hEmpty->DrawCopy();
      oneLine->Draw();
    }
    else if(iter < 8){
      hEmpty_T->DrawCopy();
      oneLine_T->Draw();
    }
    else if(iter < 12){
      hEmpty_Gen->DrawCopy();
      oneLine_Gen->Draw();
    }
    else{
      hEmpty_Gam->DrawCopy();
      oneLine_Gam->Draw();
    }

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
  for(int i = 0; i < n2Trig; ++i){
    ak2CaloTurnOns[i]->Draw("p e");
    leg[0]->AddEntry(ak2CaloTurnOns[i], ak2CaloName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 5.0");
  leg[0]->Draw("SAME");

  c1[1]->cd();
  for(int i = 0; i < n3Trig; ++i){
    ak3CaloTurnOns[i]->Draw("p e");
    leg[1]->AddEntry(ak3CaloTurnOns[i], ak3CaloName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 5.0");
  leg[1]->Draw("SAME");

  c1[2]->cd();
  for(int i = 0; i < n4Trig; ++i){
    ak4CaloTurnOns[i]->Draw("p e");
    leg[2]->AddEntry(ak4CaloTurnOns[i], ak4CaloName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 5.0");
  leg[2]->Draw("SAME");

  c1[3]->cd();
  for(int i = 0; i < n5Trig; ++i){
    ak5CaloTurnOns[i]->Draw("p e");
    leg[3]->AddEntry(ak5CaloTurnOns[i], ak5CaloName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 5.0");
  leg[3]->Draw("SAME");


  c1[4]->cd();
  for(int i = 0; i < n2Trig; ++i){
    ak2CaloTrkTurnOns[i]->Draw("p e");
    leg[4]->AddEntry(ak2CaloTrkTurnOns[i], ak2CaloTrkName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 2.4");
  leg[4]->Draw("SAME");

  c1[5]->cd();
  for(int i = 0; i < n3Trig; ++i){
    ak3CaloTrkTurnOns[i]->Draw("p e");
    leg[5]->AddEntry(ak3CaloTrkTurnOns[i], ak3CaloTrkName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 2.4");
  leg[5]->Draw("SAME");

  c1[6]->cd();
  for(int i = 0; i < n4Trig; ++i){
    ak4CaloTrkTurnOns[i]->Draw("p e");
    leg[6]->AddEntry(ak4CaloTrkTurnOns[i], ak4CaloTrkName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 2.4");
  leg[6]->Draw("SAME");

  c1[7]->cd();
  for(int i = 0; i < n5Trig; ++i){
    ak5CaloTrkTurnOns[i]->Draw("p e");
    leg[7]->AddEntry(ak5CaloTrkTurnOns[i], ak5CaloTrkName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 2.4");
  leg[7]->Draw("SAME");


  c1[8]->cd();
  for(int i = 0; i < n2Trig; ++i){
    ak2CaloGenTurnOns[i]->Draw("p e");
    leg[8]->AddEntry(ak2CaloGenTurnOns[i], ak2CaloGenName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 2.4");
  leg[8]->Draw("SAME");

  c1[9]->cd();
  for(int i = 0; i < n3Trig; ++i){
    ak3CaloGenTurnOns[i]->Draw("p e");
    leg[9]->AddEntry(ak3CaloGenTurnOns[i], ak3CaloGenName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 2.4");
  leg[9]->Draw("SAME");

  c1[10]->cd();
  for(int i = 0; i < n4Trig; ++i){
    ak4CaloGenTurnOns[i]->Draw("p e");
    leg[10]->AddEntry(ak4CaloGenTurnOns[i], ak4CaloGenName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 2.4");
  leg[10]->Draw("SAME");

  c1[11]->cd();
  for(int i = 0; i < n5Trig; ++i){
    ak5CaloGenTurnOns[i]->Draw("p e");
    leg[11]->AddEntry(ak5CaloGenTurnOns[i], ak5CaloGenName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 2.4");
  leg[11]->Draw("SAME");


  c1[12]->cd();
  for(int i = 0; i < nPhotonTrig; ++i){
    photonTurnOns[i]->Draw("p e");
    leg[12]->AddEntry(photonTurnOns[i], photonName[i], "p l");
  }
  
  label_p->DrawLatex(.75, .50, "|#eta| < 1.44");
  leg[12]->Draw("SAME");

  TFile* outFile_p = new TFile("jetTurnOnPlots_HI.root", "UPDATE");
  rateCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(rateCanv_p, Form("pdfDir/jetRateCanv"), "pdf");

  for(Int_t iter = 0; iter < nCanv; iter++){
    c1[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(c1[iter], Form("pdfDir/RECOETA5_HI%s", canvName[iter].Data()), "pdf");
  }
  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nCanv; iter++){
    delete leg[iter];
  }

  delete oneLine;

  for(Int_t iter = 0; iter < nCanv; iter++){
    delete c1[iter];
  }
  delete hEmpty;
  delete hEmpty_T;
  delete hEmpty_Gen;
  delete hEmpty_Gam;


  inFile->Close();

  return;
}
