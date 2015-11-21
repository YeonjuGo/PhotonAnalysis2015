#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>
#include <string>


void handsomeTH1(TH1* a = 0, Int_t col = 1, Float_t size = 1, Int_t markerstyle = 20)
{
  a->SetMarkerColor(col);
  a->SetMarkerSize(size);
  a->SetMarkerStyle(markerstyle);
  a->SetLineColor(col);
  a->GetYaxis()->SetTitleOffset(1.25);
  a->GetXaxis()->CenterTitle();
  a->GetYaxis()->CenterTitle();

  return;
}

void handsomeTH1N(TH1* a = 0, Int_t col = 1)
{
  handsomeTH1(a, col);
  a->Scale(1./a->GetEntries());

  return;
}


void niceTH1(TH1* uglyTH1, float max , float min, float ndivX, float ndivY, Int_t col = 1)
{
  handsomeTH1(uglyTH1, col);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->SetNdivisions(ndivX);
  uglyTH1->SetNdivisions(ndivY, "Y");

  uglyTH1->GetXaxis()->SetTitleColor(1);

  return;
}

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif")
{
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
  return;
}


void basicJetPlots(const std::string inName1, const std::string inName2, const std::string outName)
{
  TH1::SetDefaultSumw2();

  //  const std::string genLeadJtName[3] = {"genLeadJet4Pt_h", "genLeadJet4Eta_h", "genLeadJet4Phi_h"};
  //  const std::string recoLeadJtName[3] = {"recoLeadJet4Pt_h", "recoLeadJet4Eta_h", "recoLeadJet4Phi_h"};

  TCanvas* plotPanel_p = new TCanvas("plotPanel_p", "plotPanel_p", 2*300, 2*300);

  const std::string inName[2] = {inName1, inName2};
  const std::string drawOpt[4] = {"E0", "E0 SAME", "E0 SAME", "E0 SAME"};
  const Int_t color[4] = {1, kRed, kBlue, kMagenta};
  const std::string label[2] = {"twoByTwoZeroWalls", "twoByTwoZeroWallSigSub"};

  TLegend* leg_p = new TLegend(0.35, 0.65, 0.65, 0.85);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSizePixels(28);
  leg_p->SetBorderSize(0);

  for(Int_t iter = 0; iter < 2; iter++){
    TFile* inFile_p = new TFile(inName[iter].c_str(), "READ");

    TH1F* leadJetFake_h = (TH1F*)inFile_p->Get("leadJetFake_h");
    TH1F* leadJetFake2_h = (TH1F*)inFile_p->Get("leadJetFake2_h");

    niceTH1(leadJetFake_h, .99, -0.01, 505, 505, kBlue);
    niceTH1(leadJetFake2_h, .99, -0.01, 505, 505, color[iter]);
    leadJetFake2_h->SetXTitle("Scaled L1 Jet Pt");
    leadJetFake2_h->SetYTitle("Fake Rate");
    leadJetFake2_h->GetXaxis()->SetTitleOffset(1.0);

    plotPanel_p->cd();
    leadJetFake2_h->DrawCopy(drawOpt[iter].c_str());
    if(iter == 0) leadJetFake_h->DrawCopy("E0 SAME");

    if(iter == 0) leg_p->Draw("SAME");

    inFile_p->Close();
    delete inFile_p;
  }



  TH1F* dumHist_p[3];

  for(Int_t iter = 0; iter < 3; iter++){
    dumHist_p[iter] = new TH1F(Form("dum%d", iter), Form("dum%d", iter), 10, 0, 1);
    niceTH1(dumHist_p[iter], 1.20, 0.80, 505, 505, color[iter]);

    if(iter != 2) leg_p->AddEntry(dumHist_p[iter], label[iter].c_str(), "P L");
  }
  leg_p->AddEntry(dumHist_p[2], "Offline Fake", "P L");

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");
  plotPanel_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotPanel_p, "pdfDir/L1Fake2Rate", "pdf");
  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < 2; iter++){
    delete dumHist_p[iter];
  }

  delete leg_p;
  delete plotPanel_p;

  return;
}








void basicJetPlots2(const std::string inName, const std::string outName)
{
  TH1::SetDefaultSumw2();

  //  const std::string genLeadJtName[3] = {"genLeadJet4Pt_h", "genLeadJet4Eta_h", "genLeadJet4Phi_h"};
  //  const std::string recoLeadJtName[3] = {"recoLeadJet4Pt_h", "recoLeadJet4Eta_h", "recoLeadJet4Phi_h"};

  TCanvas* plotPanel_p = new TCanvas("plotPanel_p", "plotPanel_p", 2*300, 2*300);

  const std::string drawOpt[3] = {"E1", "E1 SAME", "E1 SAME"};
  const Int_t color[3] = {1, kRed, kBlue};
  const std::string label[3] = {"Single Muon", "Muon w/ Jet 60", "Muon w/ Gamma 15"};

  const std::string histName[3] = {"muRateHist_h", "muJetRateHist_h", "muGammaRateHist_h"};

  TLegend* leg_p = new TLegend(0.35, 0.65, 0.65, 0.85);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSizePixels(28);
  leg_p->SetBorderSize(0);

  for(Int_t iter = 0; iter < 3; iter++){
    TFile* inFile_p = new TFile(inName.c_str(), "READ");
    TH1F* hist_h = (TH1F*)inFile_p->Get(histName[iter].c_str());

    niceTH1(hist_h, 300, .1, 505, 505, color[iter]);
    hist_h->SetXTitle("Muon Trigger P_{T}");
    hist_h->SetYTitle("Rate (Hz)");
    hist_h->GetXaxis()->SetTitleOffset(1.0);

    plotPanel_p->cd();
    if(iter == 0) gPad->SetLogy();
    hist_h->DrawCopy(drawOpt[iter].c_str());
    if(iter == 0) leg_p->Draw("SAME");

    inFile_p->Close();
    delete inFile_p;
  }

  plotPanel_p->SetGridy();

  TH1F* dumHist_p[3];

  for(Int_t iter = 0; iter < 3; iter++){
    dumHist_p[iter] = new TH1F(Form("dum%d", iter), Form("dum%d", iter), 10, 0, 1);
    niceTH1(dumHist_p[iter], 1.20, 0.80, 505, 505, color[iter]);

    leg_p->AddEntry(dumHist_p[iter], label[iter].c_str(), "P L");
  }

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");
  plotPanel_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotPanel_p, "pdfDir/muTriggerRate", "pdf");
  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < 2; iter++){
    delete dumHist_p[iter];
  }

  delete leg_p;
  delete plotPanel_p;

  return;
}
