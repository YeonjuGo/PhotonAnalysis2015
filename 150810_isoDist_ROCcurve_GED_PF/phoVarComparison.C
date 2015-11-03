// Author Alex Barbieri
// Author modified by Yeonju 3 Aug 2015
// GEDphoton vs. old photon comparison
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TAxis.h"
#include "stdio.h"
#include "../../../HiForestAnalysis/hiForest.h"

void phoVarComparison()
{
	gStyle->SetOptStat(0);
  TChain *ged = new TChain("ggHiNtuplizer/EventTree");
  // AllQCDPhoton30 sample
  ged->Add("/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/PFphoton/CMSSW_7_5_0/src/merged_AllQCDPhoton30_standard_forest.root");

  TChain *old = new TChain("ggHiNtuplizerGED/EventTree");
  // AllQCDPhoton30 sample
  old->Add("/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/PFphoton/CMSSW_7_5_0/src/merged_AllQCDPhoton30_standard_forest.root");

  const TCut ptCut = "phoEt>30";
  const TCut etaCut = "abs(phoEta)<1.479";

  int gedEntries = ged->GetEntries();
  int oldEntries = old->GetEntries();

  TH1D *npho0 = new TH1D("npho0","# of photon;nPho",20,0,50);
  TH1D *npho1 = (TH1D*)npho0->Clone("npho1"); 
  TH1D *pt0 = new TH1D("pt0","E_{T};E_{T}^{#gamma}",20,0,150);
  TH1D *pt1 = (TH1D*)pt0->Clone("pt1");
  TH1D *eta0 = new TH1D("eta0","#eta;#eta^{#gamma}",25, -5, 5);
  TH1D *eta1 = (TH1D*)eta0->Clone("eta1");
  TH1D *phi0 = new TH1D("phi0","#phi;#phi^{#gamma}",20, -TMath::Pi(), TMath::Pi());
  TH1D *phi1 = (TH1D*)phi0->Clone("phi1");
  TH2D *map0 = new TH2D("map0","map;#eta;#phi",50, -5, 5,150, -TMath::Pi(), TMath::Pi());
  TH2D *map1 = (TH2D*)map0->Clone("map1");

  //TH2D *ptrat0 = new TH2D("ptrat0","ptrat0;p_{T}^{gen #gamma};p_{T}^{reco #gamma}",50,0,150,50,0,150);
  //TH2D *ptrat1 = (TH2D*)ptrat0->Clone("ptrat1");

  ged->Project(npho0->GetName(), "nPho");
  old->Project(npho1->GetName(),"nPho"); 
  ged->Project(pt0->GetName(), "phoEt",etaCut);
  old->Project(pt1->GetName(),"phoEt",etaCut);
  ged->Project(eta0->GetName(),"phoEta",ptCut);
  old->Project(eta1->GetName(),"phoEta",ptCut);
  ged->Project(phi0->GetName(),"phoPhi",ptCut&&etaCut);
  old->Project(phi1->GetName(),"phoPhi",ptCut&&etaCut);
  ged->Project(map0->GetName(),"phoPhi:phoEta",ptCut&&etaCut);
  old->Project(map1->GetName(),"phoPhi:phoEta",ptCut&&etaCut);

  //ged->Project(ptrat0->GetName(),"phoEt:genMatchedPt");
  //old->Project(ptrat1->GetName(),"phoEt:genMatchedPt");

  npho0->Scale(1./gedEntries);
  npho1->Scale(1./oldEntries); 
  pt0->Scale(1./gedEntries);
  pt1->Scale(1./oldEntries);
  eta0->Scale(1./gedEntries);
  eta1->Scale(1./oldEntries);
  phi0->Scale(1./gedEntries);
  phi1->Scale(1./oldEntries);
  map0->Scale(1./gedEntries);
  map1->Scale(1./oldEntries);


  TCanvas *c1 = new TCanvas();
  pt1->Draw();
  pt1->SetMarkerColor(kBlue);
  pt1->SetLineColor(kBlue);
  pt0->Draw("same");
  c1->SetLogy();

  TLegend *leg = new TLegend(0.75,0.75,0.85,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
//  leg->SetTextSize(20);
  leg->AddEntry(pt0,"GED","l");
  leg->AddEntry(pt1,"OLD","l");
  leg->Draw();

  c1->SaveAs("ged/barrel/pt.pdf");

  TCanvas *c0 = new TCanvas();
  npho1->Draw();
  npho1->SetMarkerColor(kBlue);
  npho1->SetLineColor(kBlue);
  npho0->Draw("same");
  leg->Draw();
  c0->SetLogy();
  c0->SaveAs("ged/barrel/npho_pt20.pdf");

  TCanvas *c2 = new TCanvas();
  eta1->Draw();
  eta1->SetMarkerColor(kBlue);
  eta1->SetLineColor(kBlue);
  eta0->Draw("same");
  leg->Draw();
  c2->SaveAs("ged/barrel/eta_pt20.pdf");

  TCanvas *c3 = new TCanvas();
  //phi0->SetMaximum(phi0->GetMaximum()*2);
  phi0->SetAxisRange(0.01,0.026, "Y");
  phi0->GetYaxis()->SetRangeUser(0.01,0.026);
  phi1->Draw();
  phi1->SetMarkerColor(kBlue);
  phi1->SetLineColor(kBlue);
  phi0->Draw("same");
  leg->Draw();
  c3->SaveAs("ged/barrel/phi_pt20.pdf");

  TCanvas *c4 = new TCanvas();
  map0->Draw("colz");

  TLegend *newleg = new TLegend(0.75,0.75,0.85,0.85);
  newleg->SetFillStyle(0);
  newleg->SetBorderSize(0);
  newleg->SetTextFont(43);
  newleg->SetTextSize(20);
  newleg->AddEntry(map0,"GED","");
  newleg->Draw();
  
  c4->SaveAs("ged/barrel/map_new_pt20.pdf");

  TCanvas *c5 = new TCanvas();
  map1->Draw("colz");

  TLegend *oldleg = new TLegend(0.75,0.75,0.85,0.85);
  oldleg->SetFillStyle(0);
  oldleg->SetBorderSize(0);
  oldleg->SetTextFont(43);
  oldleg->SetTextSize(20);
  oldleg->AddEntry(map1,"OLD","");
  oldleg->Draw();

  c5->SaveAs("ged/barrel/map_old_pt20.pdf");
#if 0
  TCanvas *c10 = new TCanvas();
  ptrat0->Draw("colz");
  ptrat0->ProfileX("prof0")->Draw("same");
  c10->SetLogz();
  newleg->Draw();
  c10->SaveAs("ged/barrel/corr_new_pt20.pdf");

  TCanvas *c11 = new TCanvas();
  ptrat1->Draw("colz");
  ptrat1->ProfileX("prof1")->Draw("same");
  c11->SetLogz();
  oldleg->Draw();
  c11->SaveAs("ged/barrel/corr_old_pt20.pdf");
#endif
  //isolation vars
  TH1D *sigmaIetaIeta0 = new TH1D("sigmaIetaIeta0",";#sigma_{#eta#eta}",50, 0, 0.05);
  TH1D *sigmaIetaIeta1 = (TH1D*)sigmaIetaIeta0->Clone("sigmaIetaIeta1");
  TH1D *sigmaIetaIeta20120 = new TH1D("sigmaIetaIeta20120",";#sigma_{#eta#eta}2012",50, 0, 0.05);
  TH1D *sigmaIetaIeta20121 = (TH1D*)sigmaIetaIeta0->Clone("sigmaIetaIeta20121");
  TH1D *hovere0 = new TH1D("hovere0",";h/e",50, 0, 1.2);
  TH1D *hovere1 = (TH1D*)hovere0->Clone("hovere1");
  TH1D *ecalClusterIsoR20= new TH1D("ecalClusterIsoR20",";ecalClusterIsoR2",50, -50, 50);
  TH1D *ecalClusterIsoR21= (TH1D*)ecalClusterIsoR20->Clone("ecalClusterIsoR21");
  TH1D *ecalClusterIsoR30= new TH1D("ecalClusterIsoR30",";ecalClusterIsoR3",50, -50, 50);
  TH1D *ecalClusterIsoR31= (TH1D*)ecalClusterIsoR30->Clone("ecalClusterIsoR31");
  TH1D *ecalClusterIsoR40= new TH1D("ecalClusterIsoR40",";ecalClusterIsoR4",50, -50, 50);
  TH1D *ecalClusterIsoR41= (TH1D*)ecalClusterIsoR40->Clone("ecalClusterIsoR41");
  TH1D *ecalClusterIsoR50= new TH1D("ecalClusterIsoR50",";ecalClusterIsoR5",50, -50, 50);
  TH1D *ecalClusterIsoR51= (TH1D*)ecalClusterIsoR50->Clone("ecalClusterIsoR51");
  TH1D *hcalRechitIsoR10= new TH1D("hcalRechitIsoR10",";hcalRechitIsoR1",50, -50, 50);
  TH1D *hcalRechitIsoR11= (TH1D*)hcalRechitIsoR10->Clone("hcalRechitIsoR11");
  TH1D *hcalRechitIsoR20= new TH1D("hcalRechitIsoR20",";hcalRechitIsoR2",50, -50, 50);
  TH1D *hcalRechitIsoR21= (TH1D*)hcalRechitIsoR20->Clone("hcalRechitIsoR21");
  TH1D *hcalRechitIsoR30= new TH1D("hcalRechitIsoR30",";hcalRechitIsoR3",50, -50, 50);
  TH1D *hcalRechitIsoR31= (TH1D*)hcalRechitIsoR30->Clone("hcalRechitIsoR31");
  TH1D *hcalRechitIsoR40= new TH1D("hcalRechitIsoR40",";hcalRechitIsoR4",50, -50, 50);
  TH1D *hcalRechitIsoR41= (TH1D*)hcalRechitIsoR40->Clone("hcalRechitIsoR41");
  TH1D *hcalRechitIsoR50= new TH1D("hcalRechitIsoR50",";hcalRechitIsoR5",50, -50, 50);
  TH1D *hcalRechitIsoR51= (TH1D*)hcalRechitIsoR50->Clone("hcalRechitIsoR51");
  TH1D *trackIsoR1PtCut200= new TH1D("trackIsoR1PtCut200",";trackIsoR1PtCut20",50, -20, 100);
  TH1D *trackIsoR1PtCut201= (TH1D*)trackIsoR1PtCut200->Clone("trackIsoR1PtCut201");
  TH1D *trackIsoR2PtCut200= new TH1D("trackIsoR2PtCut200",";trackIsoR2PtCut20",50, -20, 100);
  TH1D *trackIsoR2PtCut201= (TH1D*)trackIsoR2PtCut200->Clone("trackIsoR2PtCut201");
  TH1D *trackIsoR3PtCut200= new TH1D("trackIsoR3PtCut200",";trackIsoR3PtCut20",50, -20, 100);
  TH1D *trackIsoR3PtCut201= (TH1D*)trackIsoR3PtCut200->Clone("trackIsoR3PtCut201");
  TH1D *trackIsoR4PtCut200= new TH1D("trackIsoR4PtCut200",";trackIsoR4PtCut20",50, -20, 100);
  TH1D *trackIsoR4PtCut201= (TH1D*)trackIsoR4PtCut200->Clone("trackIsoR4PtCut201");
  TH1D *trackIsoR5PtCut200= new TH1D("trackIsoR5PtCut200",";trackIsoR5PtCut20",50, -20, 100);
  TH1D *trackIsoR5PtCut201= (TH1D*)trackIsoR5PtCut200->Clone("trackIsoR5PtCut201");

  ged->Project(sigmaIetaIeta0->GetName(),"phoSigmaIEtaIEta",ptCut&&etaCut);
  old->Project(sigmaIetaIeta1->GetName(),"phoSigmaIEtaIEta",ptCut&&etaCut);
  ged->Project(sigmaIetaIeta20120->GetName(),"phoSigmaIEtaIEta_2012",ptCut&&etaCut);
  old->Project(sigmaIetaIeta20121->GetName(),"phoSigmaIEtaIEta_2012",ptCut&&etaCut);
  ged->Project(hovere0->GetName(),"phoHoverE",ptCut&&etaCut);
  old->Project(hovere1->GetName(),"phoHoverE",ptCut&&etaCut);
  ged->Project(ecalClusterIsoR20->GetName(),"pho_ecalClusterIsoR2",ptCut&&etaCut);
  old->Project(ecalClusterIsoR21->GetName(),"pho_ecalClusterIsoR2",ptCut&&etaCut);
  ged->Project(ecalClusterIsoR30->GetName(),"pho_ecalClusterIsoR3",ptCut&&etaCut);
  old->Project(ecalClusterIsoR31->GetName(),"pho_ecalClusterIsoR3",ptCut&&etaCut);
  ged->Project(ecalClusterIsoR40->GetName(),"pho_ecalClusterIsoR4",ptCut&&etaCut);
  old->Project(ecalClusterIsoR41->GetName(),"pho_ecalClusterIsoR4",ptCut&&etaCut);
  ged->Project(ecalClusterIsoR50->GetName(),"pho_ecalClusterIsoR5",ptCut&&etaCut);
  old->Project(ecalClusterIsoR51->GetName(),"pho_ecalClusterIsoR5",ptCut&&etaCut);
  ged->Project(hcalRechitIsoR10->GetName(),"pho_hcalRechitIsoR1",ptCut&&etaCut);
  old->Project(hcalRechitIsoR11->GetName(),"pho_hcalRechitIsoR1",ptCut&&etaCut);
  ged->Project(hcalRechitIsoR20->GetName(),"pho_hcalRechitIsoR2",ptCut&&etaCut);
  old->Project(hcalRechitIsoR21->GetName(),"pho_hcalRechitIsoR2",ptCut&&etaCut);
  ged->Project(hcalRechitIsoR30->GetName(),"pho_hcalRechitIsoR3",ptCut&&etaCut);
  old->Project(hcalRechitIsoR31->GetName(),"pho_hcalRechitIsoR3",ptCut&&etaCut);
  ged->Project(hcalRechitIsoR40->GetName(),"pho_hcalRechitIsoR4",ptCut&&etaCut);
  old->Project(hcalRechitIsoR41->GetName(),"pho_hcalRechitIsoR4",ptCut&&etaCut);
  ged->Project(hcalRechitIsoR50->GetName(),"pho_hcalRechitIsoR5",ptCut&&etaCut);
  old->Project(hcalRechitIsoR51->GetName(),"pho_hcalRechitIsoR5",ptCut&&etaCut);
  ged->Project(trackIsoR1PtCut200->GetName(),"pho_trackIsoR1PtCut20",ptCut&&etaCut);
  old->Project(trackIsoR1PtCut201->GetName(),"pho_trackIsoR1PtCut20",ptCut&&etaCut);
  ged->Project(trackIsoR2PtCut200->GetName(),"pho_trackIsoR2PtCut20",ptCut&&etaCut);
  old->Project(trackIsoR2PtCut201->GetName(),"pho_trackIsoR2PtCut20",ptCut&&etaCut);
  ged->Project(trackIsoR3PtCut200->GetName(),"pho_trackIsoR3PtCut20",ptCut&&etaCut);
  old->Project(trackIsoR3PtCut201->GetName(),"pho_trackIsoR3PtCut20",ptCut&&etaCut);
  ged->Project(trackIsoR4PtCut200->GetName(),"pho_trackIsoR4PtCut20",ptCut&&etaCut);
  old->Project(trackIsoR4PtCut201->GetName(),"pho_trackIsoR4PtCut20",ptCut&&etaCut);
  ged->Project(trackIsoR5PtCut200->GetName(),"pho_trackIsoR5PtCut20",ptCut&&etaCut);
  old->Project(trackIsoR5PtCut201->GetName(),"pho_trackIsoR5PtCut20",ptCut&&etaCut);
 #if 0 
  cc40->Scale(1./gedEntries);
  cc41->Scale(1./oldEntries);
#endif

  TCanvas *c6 = new TCanvas();
  sigmaIetaIeta0->Draw();
  sigmaIetaIeta1->SetMarkerColor(kBlue);
  sigmaIetaIeta1->SetLineColor(kBlue);
  sigmaIetaIeta1->Draw("same");
  c6->SetLogy();
  leg->Draw();
  c6->SaveAs("ged/barrel/sigmaIetaIeta_pt20.pdf");

  TCanvas *c6a = new TCanvas();
  sigmaIetaIeta20120->Draw();
  sigmaIetaIeta20121->SetMarkerColor(kBlue);
  sigmaIetaIeta20121->SetLineColor(kBlue);
  sigmaIetaIeta20121->Draw("same");
  c6a->SetLogy();
  leg->Draw();
  c6a->SaveAs("ged/barrel/sigmaIetaIeta2012_pt20.pdf");

  TCanvas *c7 = new TCanvas();
  hovere0->Draw();
  hovere1->SetMarkerColor(kBlue);
  hovere1->SetLineColor(kBlue);
  hovere1->Draw("same");
  c7->SetLogy();
  leg->Draw();
  c7->SaveAs("ged/barrel/hovere_pt20.pdf");
#if 0
  TCanvas *c8 = new TCanvas();
  ecalClusterIsoR20->Draw();
  ecalClusterIsoR21->SetMarkerColor(kBlue);
  ecalClusterIsoR21->SetLineColor(kBlue);
  ecalClusterIsoR21->Draw("same");
  c8->SetLogy();
  leg->Draw();
  c8->SaveAs("ged/barrel/ecalClusterIsoR2_pt20.pdf");

  TCanvas *c8a = new TCanvas();
  ecalClusterIsoR30->Draw();
  ecalClusterIsoR31->SetMarkerColor(kBlue);
  ecalClusterIsoR31->SetLineColor(kBlue);
  ecalClusterIsoR31->Draw("same");
  c8a->SetLogy();
  leg->Draw();
  c8a->SaveAs("ged/barrel/ecalClusterIsoR3_pt20.pdf");

  TCanvas *c8b = new TCanvas();
  ecalClusterIsoR40->Draw();
  ecalClusterIsoR41->SetMarkerColor(kBlue);
  ecalClusterIsoR41->SetLineColor(kBlue);
  ecalClusterIsoR41->Draw("same");
  c8b->SetLogy();
  leg->Draw();
  c8b->SaveAs("ged/barrel/ecalClusterIsoR4_pt20.pdf");

  TCanvas *c8c = new TCanvas();
  ecalClusterIsoR50->Draw();
  ecalClusterIsoR51->SetMarkerColor(kBlue);
  ecalClusterIsoR51->SetLineColor(kBlue);
  ecalClusterIsoR51->Draw("same");
  c8c->SetLogy();
  leg->Draw();
  c8c->SaveAs("ged/barrel/ecalClusterIsoR5_pt20.pdf");

  TCanvas *c9 = new TCanvas();
  hcalRechitIsoR10->Draw();
  hcalRechitIsoR11->SetMarkerColor(kBlue);
  hcalRechitIsoR11->SetLineColor(kBlue);
  hcalRechitIsoR11->Draw("same");
  c9->SetLogy();
  leg->Draw();
  c9->SaveAs("ged/barrel/hcalRechitIsoR1_pt20.pdf");

  TCanvas *c9a = new TCanvas();
  hcalRechitIsoR20->Draw();
  hcalRechitIsoR21->SetMarkerColor(kBlue);
  hcalRechitIsoR21->SetLineColor(kBlue);
  hcalRechitIsoR21->Draw("same");
  c9a->SetLogy();
  leg->Draw();
  c9a->SaveAs("ged/barrel/hcalRechitIsoR2_pt20.pdf");

  TCanvas *c9b = new TCanvas();
  hcalRechitIsoR30->Draw();
  hcalRechitIsoR31->SetMarkerColor(kBlue);
  hcalRechitIsoR31->SetLineColor(kBlue);
  hcalRechitIsoR31->Draw("same");
  c9b->SetLogy();
  leg->Draw();
  c9b->SaveAs("ged/barrel/hcalRechitIsoR3_pt20.pdf");

  TCanvas *c9c = new TCanvas();
  hcalRechitIsoR40->Draw();
  hcalRechitIsoR41->SetMarkerColor(kBlue);
  hcalRechitIsoR41->SetLineColor(kBlue);
  hcalRechitIsoR41->Draw("same");
  c9c->SetLogy();
  leg->Draw();
  c9c->SaveAs("ged/barrel/hcalRechitIsoR4_pt20.pdf");

  TCanvas *c9d = new TCanvas();
  hcalRechitIsoR50->Draw();
  hcalRechitIsoR51->SetMarkerColor(kBlue);
  hcalRechitIsoR51->SetLineColor(kBlue);
  hcalRechitIsoR51->Draw("same");
  c9d->SetLogy();
  leg->Draw();
  c9d->SaveAs("ged/barrel/hcalRechitIsoR5_pt20.pdf");

  TCanvas *c10 = new TCanvas();
  trackIsoR1PtCut200->Draw();
  trackIsoR1PtCut201->SetMarkerColor(kBlue);
  trackIsoR1PtCut201->SetLineColor(kBlue);
  trackIsoR1PtCut201->Draw("same");
  c10->SetLogy();
  leg->Draw();
  c10->SaveAs("ged/barrel/trackIsoR1PtCut20_pt20.pdf");

  TCanvas *c10a = new TCanvas();
  trackIsoR2PtCut200->Draw();
  trackIsoR2PtCut201->SetMarkerColor(kBlue);
  trackIsoR2PtCut201->SetLineColor(kBlue);
  trackIsoR2PtCut201->Draw("same");
  c10a->SetLogy();
  leg->Draw();
  c10a->SaveAs("ged/barrel/trackIsoR2PtCut20_pt20.pdf");

  TCanvas *c10b = new TCanvas();
  trackIsoR3PtCut200->Draw();
  trackIsoR3PtCut201->SetMarkerColor(kBlue);
  trackIsoR3PtCut201->SetLineColor(kBlue);
  trackIsoR3PtCut201->Draw("same");
  c10b->SetLogy();
  leg->Draw();
  c10b->SaveAs("ged/barrel/trackIsoR3PtCut20_pt20.pdf");

  TCanvas *c10c = new TCanvas();
  trackIsoR4PtCut200->Draw();
  trackIsoR4PtCut201->SetMarkerColor(kBlue);
  trackIsoR4PtCut201->SetLineColor(kBlue);
  trackIsoR4PtCut201->Draw("same");
  c10c->SetLogy();
  leg->Draw();
  c10c->SaveAs("ged/barrel/trackIsoR4PtCut20_pt20.pdf");

  TCanvas *c10d = new TCanvas();
  trackIsoR5PtCut200->Draw();
  trackIsoR5PtCut201->SetMarkerColor(kBlue);
  trackIsoR5PtCut201->SetLineColor(kBlue);
  trackIsoR5PtCut201->Draw("same");
  c10d->SetLogy();
  leg->Draw();
  c10d->SaveAs("ged/barrel/trackIsoR5PtCut20_pt20.pdf");
#endif
}
  
