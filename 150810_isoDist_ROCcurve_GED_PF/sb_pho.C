// Author Alex Barbieri
// Author modified by Yeonju 3 Aug 2015
// GEDphoton vs. sig photon comparison
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

void sb_pho(bool isGED=false, isGenMatch=true)
{
  gStyle->SetOptStat(0);
  gStyle->SetHistLineWidth(2);

  TString st="";
  if(isGED) st="GED";
  TChain *bkg = new TChain(Form("ggHiNtuplizer%s/EventTree",st.Data()));
  // EmEnrichedDijet30 sample
  bkg->Add("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_EmEnrichedDijet30_standard_forest.root");

  TChain *sig = new TChain(Form("ggHiNtuplizer%s/EventTree",st.Data()));
  // AllQCDPhoton30 sample
  sig->Add("/afs/cern.ch/work/y/ygo/public/PFphoton/merged_AllQCDPhoton30_standard_forest_2nd.root");

  int ptval = 30;
  const TCut ptCut = Form("phoEt>%d",ptval);
  const TCut etaCut = "abs(phoEta)<1.479";

  int bkgEntries = bkg->GetEntries();
  int sigEntries = sig->GetEntries();

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

  bkg->Project(npho0->GetName(), "nPho");
  sig->Project(npho1->GetName(),"nPho"); 
  bkg->Project(pt0->GetName(), "phoEt",etaCut);
  sig->Project(pt1->GetName(),"phoEt",etaCut);
  bkg->Project(eta0->GetName(),"phoEta",ptCut);
  sig->Project(eta1->GetName(),"phoEta",ptCut);
  bkg->Project(phi0->GetName(),"phoPhi",ptCut&&etaCut);
  sig->Project(phi1->GetName(),"phoPhi",ptCut&&etaCut);
  bkg->Project(map0->GetName(),"phoPhi:phoEta",ptCut&&etaCut);
  sig->Project(map1->GetName(),"phoPhi:phoEta",ptCut&&etaCut);

  //bkg->Project(ptrat0->GetName(),"phoEt:genMatchedPt");
  //sig->Project(ptrat1->GetName(),"phoEt:genMatchedPt");

  npho0->Scale(1./bkgEntries);
  npho1->Scale(1./sigEntries); 
  pt0->Scale(1./bkgEntries);
  pt1->Scale(1./sigEntries);
  eta0->Scale(1./bkgEntries);
  eta1->Scale(1./sigEntries);
  phi0->Scale(1./bkgEntries);
  phi1->Scale(1./sigEntries);
  map0->Scale(1./bkgEntries);
  map1->Scale(1./sigEntries);


  TCanvas *c1 = new TCanvas();
  pt0->Draw();
  pt1->SetMarkerColor(kBlue);
  pt1->SetLineColor(kBlue);
  pt1->Draw("same");
  c1->SetLogy();

  TLegend *leg = new TLegend(0.60,0.75,0.85,0.89);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
//  leg->SetTextSize(20);
  leg->AddEntry(pt0,"EmEnrichedDijet30","l");
  leg->AddEntry(pt1,"AllQCDPhoton30","l");
  leg->Draw();

  c1->SaveAs(Form("bkg/barrel/%s_pt.pdf",st.Data()));

  TCanvas *c0 = new TCanvas();
  npho0->Draw();
  npho1->SetMarkerColor(kBlue);
  npho1->SetLineColor(kBlue);
  npho1->Draw("same");
  leg->Draw();
  c0->SetLogy();
  c0->SaveAs(Form("bkg/barrel/%s_npho.pdf",st.Data()));

  TCanvas *c2 = new TCanvas();
  eta1->Draw();
  eta1->SetMarkerColor(kBlue);
  eta1->SetLineColor(kBlue);
  eta0->Draw("same");
  leg->Draw();
  c2->SaveAs(Form("bkg/barrel/%s_eta_pt%d.pdf",st.Data(),ptval));

  TCanvas *c3 = new TCanvas();
  //phi0->SetMaximum(phi0->GetMaximum()*2);
  phi0->SetAxisRange(0.0,0.026, "Y");
  phi1->SetAxisRange(0.0,0.026, "Y");
 // phi0->GetYaxis()->SetRangeUser(0.01,0.026);
  phi1->Draw();
  phi1->SetMarkerColor(kBlue);
  phi1->SetLineColor(kBlue);
  phi0->Draw("same");
  leg->Draw();
  c3->SaveAs(Form("bkg/barrel/%s_phi_pt%d.pdf",st.Data(),ptval));

#if 0
  TCanvas *c4 = new TCanvas();
  map0->Draw("colz");

  TLegend *newleg = new TLegend(0.65,0.75,0.85,0.85);
  newleg->SetFillStyle(0);
  newleg->SetBorderSize(0);
  newleg->SetTextFont(43);
  newleg->SetTextSize(20);
  newleg->AddEntry(map0,"EmEnrichedDijet30","");
  newleg->Draw();
  
  c4->SaveAs(Form("bkg/barrel/%s_map_bkg_pt%d.pdf",st.Data(),ptval));

  TCanvas *c5 = new TCanvas();
  map1->Draw("colz");

  TLegend *sigleg = new TLegend(0.65,0.75,0.85,0.85);
  sigleg->SetFillStyle(0);
  sigleg->SetBorderSize(0);
  sigleg->SetTextFont(43);
  sigleg->SetTextSize(20);
  sigleg->AddEntry(map1,"AllQCDPhoton30","");
  sigleg->Draw();

  c5->SaveAs(Form("bkg/barrel/%s_map_sig_pt%d.pdf",st.Data(),ptval));

  TCanvas *c10 = new TCanvas();
  ptrat0->Draw("colz");
  ptrat0->ProfileX("prof0")->Draw("same");
  c10->SetLogz();
  newleg->Draw();
  c10->SaveAs("bkg/barrel/corr_new_pt20.pdf");

  TCanvas *c11 = new TCanvas();
  ptrat1->Draw("colz");
  ptrat1->ProfileX("prof1")->Draw("same");
  c11->SetLogz();
  sigleg->Draw();
  c11->SaveAs("bkg/barrel/corr_sig_pt20.pdf");
#endif
  //isolation vars
  TH1D *sigmaIetaIeta0 = new TH1D("sigmaIetaIeta0","#sigma_{I#etaI#eta};#sigma_{#eta#eta}",50, 0, 0.05);
  TH1D *sigmaIetaIeta1 = (TH1D*)sigmaIetaIeta0->Clone("sigmaIetaIeta1");
  TH1D *sigmaIetaIeta20120 = new TH1D("sigmaIetaIeta20120",";#sigma_{#eta#eta}2012",50, 0, 0.05);
  TH1D *sigmaIetaIeta20121 = (TH1D*)sigmaIetaIeta0->Clone("sigmaIetaIeta20121");
  TH1D *hovere0 = new TH1D("hovere0","h/e;h/e",50, 0, 1.2);
  TH1D *hovere1 = (TH1D*)hovere0->Clone("hovere1");
  TH1D *ecalClusterIsoR20= new TH1D("ecalClusterIsoR20",";ecalClusterIsoR2",50, -40, 90);
  TH1D *ecalClusterIsoR21= (TH1D*)ecalClusterIsoR20->Clone("ecalClusterIsoR21");
  TH1D *ecalClusterIsoR30= new TH1D("ecalClusterIsoR30",";ecalClusterIsoR3",50, -40, 90);
  TH1D *ecalClusterIsoR31= (TH1D*)ecalClusterIsoR30->Clone("ecalClusterIsoR31");
  TH1D *ecalClusterIsoR40= new TH1D("ecalClusterIsoR40",";ecalClusterIsoR4",50, -40, 90);
  TH1D *ecalClusterIsoR41= (TH1D*)ecalClusterIsoR40->Clone("ecalClusterIsoR41");
  TH1D *ecalClusterIsoR50= new TH1D("ecalClusterIsoR50",";ecalClusterIsoR5",50, -40, 90);
  TH1D *ecalClusterIsoR51= (TH1D*)ecalClusterIsoR50->Clone("ecalClusterIsoR51");
  TH1D *hcalRechitIsoR10= new TH1D("hcalRechitIsoR10",";hcalRechitIsoR1",50, -40, 90);
  TH1D *hcalRechitIsoR11= (TH1D*)hcalRechitIsoR10->Clone("hcalRechitIsoR11");
  TH1D *hcalRechitIsoR20= new TH1D("hcalRechitIsoR20",";hcalRechitIsoR2",50, -40, 90);
  TH1D *hcalRechitIsoR21= (TH1D*)hcalRechitIsoR20->Clone("hcalRechitIsoR21");
  TH1D *hcalRechitIsoR30= new TH1D("hcalRechitIsoR30",";hcalRechitIsoR3",50, -40, 90);
  TH1D *hcalRechitIsoR31= (TH1D*)hcalRechitIsoR30->Clone("hcalRechitIsoR31");
  TH1D *hcalRechitIsoR40= new TH1D("hcalRechitIsoR40",";hcalRechitIsoR4",50, -40, 90);
  TH1D *hcalRechitIsoR41= (TH1D*)hcalRechitIsoR40->Clone("hcalRechitIsoR41");
  TH1D *hcalRechitIsoR50= new TH1D("hcalRechitIsoR50",";hcalRechitIsoR5",50, -40, 90);
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

  bkg->Project(sigmaIetaIeta0->GetName(),"phoSigmaIEtaIEta",ptCut&&etaCut);
  sig->Project(sigmaIetaIeta1->GetName(),"phoSigmaIEtaIEta",ptCut&&etaCut);
  bkg->Project(sigmaIetaIeta20120->GetName(),"phoSigmaIEtaIEta_2012",ptCut&&etaCut);
  sig->Project(sigmaIetaIeta20121->GetName(),"phoSigmaIEtaIEta_2012",ptCut&&etaCut);
  bkg->Project(hovere0->GetName(),"phoHoverE",ptCut&&etaCut);
  sig->Project(hovere1->GetName(),"phoHoverE",ptCut&&etaCut);
  bkg->Project(ecalClusterIsoR20->GetName(),"pho_ecalClusterIsoR2",ptCut&&etaCut);
  sig->Project(ecalClusterIsoR21->GetName(),"pho_ecalClusterIsoR2",ptCut&&etaCut);
  bkg->Project(ecalClusterIsoR30->GetName(),"pho_ecalClusterIsoR3",ptCut&&etaCut);
  sig->Project(ecalClusterIsoR31->GetName(),"pho_ecalClusterIsoR3",ptCut&&etaCut);
  bkg->Project(ecalClusterIsoR40->GetName(),"pho_ecalClusterIsoR4",ptCut&&etaCut);
  sig->Project(ecalClusterIsoR41->GetName(),"pho_ecalClusterIsoR4",ptCut&&etaCut);
  bkg->Project(ecalClusterIsoR50->GetName(),"pho_ecalClusterIsoR5",ptCut&&etaCut);
  sig->Project(ecalClusterIsoR51->GetName(),"pho_ecalClusterIsoR5",ptCut&&etaCut);
  bkg->Project(hcalRechitIsoR10->GetName(),"pho_hcalRechitIsoR1",ptCut&&etaCut);
  sig->Project(hcalRechitIsoR11->GetName(),"pho_hcalRechitIsoR1",ptCut&&etaCut);
  bkg->Project(hcalRechitIsoR20->GetName(),"pho_hcalRechitIsoR2",ptCut&&etaCut);
  sig->Project(hcalRechitIsoR21->GetName(),"pho_hcalRechitIsoR2",ptCut&&etaCut);
  bkg->Project(hcalRechitIsoR30->GetName(),"pho_hcalRechitIsoR3",ptCut&&etaCut);
  sig->Project(hcalRechitIsoR31->GetName(),"pho_hcalRechitIsoR3",ptCut&&etaCut);
  bkg->Project(hcalRechitIsoR40->GetName(),"pho_hcalRechitIsoR4",ptCut&&etaCut);
  sig->Project(hcalRechitIsoR41->GetName(),"pho_hcalRechitIsoR4",ptCut&&etaCut);
  bkg->Project(hcalRechitIsoR50->GetName(),"pho_hcalRechitIsoR5",ptCut&&etaCut);
  sig->Project(hcalRechitIsoR51->GetName(),"pho_hcalRechitIsoR5",ptCut&&etaCut);
  bkg->Project(trackIsoR1PtCut200->GetName(),"pho_trackIsoR1PtCut20",ptCut&&etaCut);
  sig->Project(trackIsoR1PtCut201->GetName(),"pho_trackIsoR1PtCut20",ptCut&&etaCut);
  bkg->Project(trackIsoR2PtCut200->GetName(),"pho_trackIsoR2PtCut20",ptCut&&etaCut);
  sig->Project(trackIsoR2PtCut201->GetName(),"pho_trackIsoR2PtCut20",ptCut&&etaCut);
  bkg->Project(trackIsoR3PtCut200->GetName(),"pho_trackIsoR3PtCut20",ptCut&&etaCut);
  sig->Project(trackIsoR3PtCut201->GetName(),"pho_trackIsoR3PtCut20",ptCut&&etaCut);
  bkg->Project(trackIsoR4PtCut200->GetName(),"pho_trackIsoR4PtCut20",ptCut&&etaCut);
  sig->Project(trackIsoR4PtCut201->GetName(),"pho_trackIsoR4PtCut20",ptCut&&etaCut);
  bkg->Project(trackIsoR5PtCut200->GetName(),"pho_trackIsoR5PtCut20",ptCut&&etaCut);
  sig->Project(trackIsoR5PtCut201->GetName(),"pho_trackIsoR5PtCut20",ptCut&&etaCut);

  sigmaIetaIeta0->Scale(1./bkgEntries);
  sigmaIetaIeta1->Scale(1./sigEntries);
  sigmaIetaIeta20120->Scale(1./bkgEntries);
  sigmaIetaIeta20121->Scale(1./sigEntries);
  hovere0->Scale(1./bkgEntries);
  hovere1->Scale(1./sigEntries);
  ecalClusterIsoR20->Scale(1./bkgEntries);
  ecalClusterIsoR21->Scale(1./sigEntries);
  ecalClusterIsoR30->Scale(1./bkgEntries);
  ecalClusterIsoR31->Scale(1./sigEntries);
  ecalClusterIsoR40->Scale(1./bkgEntries);
  ecalClusterIsoR41->Scale(1./sigEntries);
  ecalClusterIsoR50->Scale(1./bkgEntries);
  ecalClusterIsoR51->Scale(1./sigEntries);
  hcalRechitIsoR10->Scale(1./bkgEntries);
  hcalRechitIsoR11->Scale(1./sigEntries);
  hcalRechitIsoR20->Scale(1./bkgEntries);
  hcalRechitIsoR21->Scale(1./sigEntries);
  hcalRechitIsoR30->Scale(1./bkgEntries);
  hcalRechitIsoR31->Scale(1./sigEntries);
  hcalRechitIsoR40->Scale(1./bkgEntries);
  hcalRechitIsoR41->Scale(1./sigEntries);
  hcalRechitIsoR50->Scale(1./bkgEntries);
  hcalRechitIsoR51->Scale(1./sigEntries);
  trackIsoR1PtCut200->Scale(1./bkgEntries);
  trackIsoR1PtCut201->Scale(1./sigEntries);
  trackIsoR2PtCut200->Scale(1./bkgEntries);
  trackIsoR2PtCut201->Scale(1./sigEntries);
  trackIsoR3PtCut200->Scale(1./bkgEntries);
  trackIsoR3PtCut201->Scale(1./sigEntries);
  trackIsoR4PtCut200->Scale(1./bkgEntries);
  trackIsoR4PtCut201->Scale(1./sigEntries);
  trackIsoR5PtCut200->Scale(1./bkgEntries);
  trackIsoR5PtCut201->Scale(1./sigEntries);



  TCanvas *c6 = new TCanvas();
  sigmaIetaIeta1->Draw();
  sigmaIetaIeta1->SetMarkerColor(kBlue);
  sigmaIetaIeta1->SetLineColor(kBlue);
  sigmaIetaIeta0->Draw("same");
  c6->SetLogy();
  leg->Draw();
  c6->SaveAs(Form("bkg/barrel/%s_sigmaIetaIeta_pt%d.pdf",st.Data(),ptval));
#if 0 
  TCanvas *c6a = new TCanvas();
  sigmaIetaIeta20121->Draw();
  sigmaIetaIeta20121->SetMarkerColor(kBlue);
  sigmaIetaIeta20121->SetLineColor(kBlue);
  sigmaIetaIeta20120->Draw("same");
  c6a->SetLogy();
  leg->Draw();
  c6a->SaveAs(Form("bkg/barrel/%s_sigmaIetaIeta2012_pt%d.pdf",st.Data(),ptval));
#endif
  TCanvas *c7 = new TCanvas();
  hovere1->Draw();
  hovere1->SetMarkerColor(kBlue);
  hovere1->SetLineColor(kBlue);
  hovere0->Draw("same");
  c7->SetLogy();
  leg->Draw();
  c7->SaveAs(Form("bkg/barrel/%s_hovere_pt%d.pdf",st.Data(),ptval));
#if 1
  TCanvas *c8 = new TCanvas();
  ecalClusterIsoR21->Draw();
  ecalClusterIsoR21->SetMarkerColor(kBlue);
  ecalClusterIsoR21->SetLineColor(kBlue);
  ecalClusterIsoR20->Draw("same");
  c8->SetLogy();
  leg->Draw();
  c8->SaveAs(Form("bkg/barrel/%s_ecalClusterIsoR2_pt%d.pdf",st.Data(),ptval));

  TCanvas *c8a = new TCanvas();
  ecalClusterIsoR31->Draw();
  ecalClusterIsoR31->SetMarkerColor(kBlue);
  ecalClusterIsoR31->SetLineColor(kBlue);
  ecalClusterIsoR30->Draw("same");
  c8a->SetLogy();
  leg->Draw();
  c8a->SaveAs(Form("bkg/barrel/%s_ecalClusterIsoR3_pt%d.pdf",st.Data(),ptval));

  TCanvas *c8b = new TCanvas();
  ecalClusterIsoR41->Draw();
  ecalClusterIsoR41->SetMarkerColor(kBlue);
  ecalClusterIsoR41->SetLineColor(kBlue);
  ecalClusterIsoR40->Draw("same");
  c8b->SetLogy();
  leg->Draw();
  c8b->SaveAs(Form("bkg/barrel/%s_ecalClusterIsoR4_pt%d.pdf",st.Data(),ptval));

  TCanvas *c8c = new TCanvas();
  ecalClusterIsoR51->Draw();
  ecalClusterIsoR51->SetMarkerColor(kBlue);
  ecalClusterIsoR51->SetLineColor(kBlue);
  ecalClusterIsoR50->Draw("same");
  c8c->SetLogy();
  leg->Draw();
  c8c->SaveAs(Form("bkg/barrel/%s_ecalClusterIsoR5_pt%d.pdf",st.Data(),ptval));

  TCanvas *c9 = new TCanvas();
  hcalRechitIsoR11->Draw();
  hcalRechitIsoR11->SetMarkerColor(kBlue);
  hcalRechitIsoR11->SetLineColor(kBlue);
  hcalRechitIsoR10->Draw("same");
  c9->SetLogy();
  leg->Draw();
  c9->SaveAs(Form("bkg/barrel/%s_hcalRechitIsoR1_pt%d.pdf",st.Data(),ptval));

  TCanvas *c9a = new TCanvas();
  hcalRechitIsoR21->Draw();
  hcalRechitIsoR21->SetMarkerColor(kBlue);
  hcalRechitIsoR21->SetLineColor(kBlue);
  hcalRechitIsoR20->Draw("same");
  c9a->SetLogy();
  leg->Draw();
  c9a->SaveAs(Form("bkg/barrel/%s_hcalRechitIsoR2_pt%d.pdf",st.Data(),ptval));

  TCanvas *c9b = new TCanvas();
  hcalRechitIsoR31->Draw();
  hcalRechitIsoR31->SetMarkerColor(kBlue);
  hcalRechitIsoR31->SetLineColor(kBlue);
  hcalRechitIsoR30->Draw("same");
  c9b->SetLogy();
  leg->Draw();
  c9b->SaveAs(Form("bkg/barrel/%s_hcalRechitIsoR3_pt%d.pdf",st.Data(),ptval));

  TCanvas *c9c = new TCanvas();
  hcalRechitIsoR41->Draw();
  hcalRechitIsoR41->SetMarkerColor(kBlue);
  hcalRechitIsoR41->SetLineColor(kBlue);
  hcalRechitIsoR40->Draw("same");
  c9c->SetLogy();
  leg->Draw();
  c9c->SaveAs(Form("bkg/barrel/%s_hcalRechitIsoR4_pt%d.pdf",st.Data(),ptval));

  TCanvas *c9d = new TCanvas();
  hcalRechitIsoR51->Draw();
  hcalRechitIsoR51->SetMarkerColor(kBlue);
  hcalRechitIsoR51->SetLineColor(kBlue);
  hcalRechitIsoR50->Draw("same");
  c9d->SetLogy();
  leg->Draw();
  c9d->SaveAs(Form("bkg/barrel/%s_hcalRechitIsoR5_pt%d.pdf",st.Data(),ptval));

  TCanvas *c10 = new TCanvas();
  trackIsoR1PtCut201->Draw();
  trackIsoR1PtCut201->SetMarkerColor(kBlue);
  trackIsoR1PtCut201->SetLineColor(kBlue);
  trackIsoR1PtCut200->Draw("same");
  c10->SetLogy();
  leg->Draw();
  c10->SaveAs(Form("bkg/barrel/%s_trackIsoR1PtCut20_pt%d.pdf",st.Data(),ptval));

  TCanvas *c10a = new TCanvas();
  trackIsoR2PtCut201->Draw();
  trackIsoR2PtCut201->SetMarkerColor(kBlue);
  trackIsoR2PtCut201->SetLineColor(kBlue);
  trackIsoR2PtCut200->Draw("same");
  c10a->SetLogy();
  leg->Draw();
  c10a->SaveAs(Form("bkg/barrel/%s_trackIsoR2PtCut20_pt%d.pdf",st.Data(),ptval));

  TCanvas *c10b = new TCanvas();
  trackIsoR3PtCut201->Draw();
  trackIsoR3PtCut201->SetMarkerColor(kBlue);
  trackIsoR3PtCut201->SetLineColor(kBlue);
  trackIsoR3PtCut200->Draw("same");
  c10b->SetLogy();
  leg->Draw();
  c10b->SaveAs(Form("bkg/barrel/%s_trackIsoR3PtCut20_pt%d.pdf",st.Data(),ptval));

  TCanvas *c10c = new TCanvas();
  trackIsoR4PtCut201->Draw();
  trackIsoR4PtCut201->SetMarkerColor(kBlue);
  trackIsoR4PtCut201->SetLineColor(kBlue);
  trackIsoR4PtCut200->Draw("same");
  c10c->SetLogy();
  leg->Draw();
  c10c->SaveAs(Form("bkg/barrel/%s_trackIsoR4PtCut20_pt%d.pdf",st.Data(),ptval));

  TCanvas *c10d = new TCanvas();
  trackIsoR5PtCut201->Draw();
  trackIsoR5PtCut201->SetMarkerColor(kBlue);
  trackIsoR5PtCut201->SetLineColor(kBlue);
  trackIsoR5PtCut200->Draw("same");
  c10d->SetLogy();
  leg->Draw();
  c10d->SaveAs(Form("bkg/barrel/%s_trackIsoR5PtCut20_pt%d.pdf",st.Data(),ptval));
#endif
}
  
