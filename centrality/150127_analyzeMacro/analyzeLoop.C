#include "TNtuple.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;

#define MAXHITS 1000000

struct MyRecHit{

   int n;

   float e[MAXHITS];
   float et[MAXHITS];
   float eta[MAXHITS];
   float phi[MAXHITS];
   bool isjet[MAXHITS];

   float jtpt;
   float jteta;
   float jtphi;

   int depth;

};

struct MyBkg{
   int n;
   float rho[50];
   float sigma[50];
};


int analyzeLoop(const char* infile = "test.root", const char* outfile = "output2.root"){

   cout<<"a"<<endl;

   const char* trig = "L1Tech_BSC_minBias_threshold1.v0";
   double towerCut = 3;

   bool MC = true;

   TFile * inf = new TFile(infile);
   TFile* outf = new TFile(outfile,"recreate");

   TNtuple* nt = new TNtuple("evNT","","ETmyTowers:NmyTowers:EThfRecHit:NhfRecHit:NgenParts:ETgenParts");

   cout<<"a"<<endl;


   TH1D* ha[100];
   TH1D* hb[100];
   TH1D* hc[100];
   TH1D* he[100];
   TH1D* hn[100];

   TH2D* h2o[100];
   TH2D* h2a[100];
   TH2D* h2b[100];
   TH2D* h2c[100];
   TH2D* h2d[100];

   for(int i = 0; i < 20; ++i){
      ha[i] = new TH1D(Form("ha%d",i),"",1000,0,1000);

      hb[i] = new TH1D(Form("hb%d",i),"",1000,0,1000);
      hc[i] = new TH1D(Form("hc%d",i),"",1000,0,1000);

      he[i] = new TH1D(Form("he%d",i),"",1000,0,4000);
      hn[i] = new TH1D(Form("hn%d",i),"",1000,0,10000);

      h2o[i] = new TH2D(Form("h2o%d",i),"",420,0,420,600,0,6000);
      h2a[i] = new TH2D(Form("h2a%d",i),"",420,0,420,900,0,9000);
      h2b[i] = new TH2D(Form("h2b%d",i),"",900,0,9000,600,0,6000);
      h2c[i] = new TH2D(Form("h2c%d",i),"",420,0,420,370,0,3700);
      h2d[i] = new TH2D(Form("h2d%d",i),"",370,0,3700,600,0,6000);

   }

   ha[1]->SetLineColor(2);
   ha[1]->SetMarkerColor(2);

   TTree* t1 = (TTree*)inf->Get("hltanalysis/HltTree");
   TTree* t2 = (TTree*)inf->Get("rechits/hf");
   TTree* t3 = (TTree*)inf->Get("rechits/hbhe");
   TTree* t4 = (TTree*)inf->Get("rechits/ee");
   TTree* t5 = (TTree*)inf->Get("rechits/eb");
   TTree* t6 = (TTree*)inf->Get("rechits/bkg");
   TTree* t7 = (TTree*)inf->Get("rechits/tower");
   TTree* t8;
   if(MC){
      t8 = (TTree*)inf->Get("mc/hi");
      t1->AddFriend(t8);
   }

   cout<<"a"<<endl;
   t1->AddFriend(t2);
   t1->AddFriend(t3);
   t1->AddFriend(t4);
   t1->AddFriend(t5);

   MyRecHit hbheRecHit;
   MyRecHit hfRecHit;
   MyRecHit ebRecHit;
   MyRecHit eeRecHit;
   MyRecHit myBC;
   MyRecHit myTowers;
   MyRecHit genParts;
   MyBkg bkg;

   t3->SetBranchAddress("e",hbheRecHit.e);
   t3->SetBranchAddress("et",hbheRecHit.et);
   t3->SetBranchAddress("eta",hbheRecHit.eta);
   t3->SetBranchAddress("phi",hbheRecHit.phi);
   t3->SetBranchAddress("n",&hbheRecHit.n);

   t2->SetBranchAddress("e",hfRecHit.e);
   t2->SetBranchAddress("et",hfRecHit.et);
   t2->SetBranchAddress("eta",hfRecHit.eta);
   t2->SetBranchAddress("phi",hfRecHit.phi);
   t2->SetBranchAddress("n",&hfRecHit.n);

   t5->SetBranchAddress("e",ebRecHit.e);
   t5->SetBranchAddress("et",ebRecHit.et);
   t5->SetBranchAddress("eta",ebRecHit.eta);
   t5->SetBranchAddress("phi",ebRecHit.phi);
   t5->SetBranchAddress("n",&ebRecHit.n);

   t4->SetBranchAddress("e",eeRecHit.e);
   t4->SetBranchAddress("et",eeRecHit.et);
   t4->SetBranchAddress("eta",eeRecHit.eta);
   t4->SetBranchAddress("phi",eeRecHit.phi);
   t4->SetBranchAddress("n",&eeRecHit.n);

   t7->SetBranchAddress("e",myTowers.e);
   t7->SetBranchAddress("et",myTowers.et);
   t7->SetBranchAddress("eta",myTowers.eta);
   t7->SetBranchAddress("phi",myTowers.phi);
   t7->SetBranchAddress("n",&myTowers.n);

   t8->SetBranchAddress("pt",genParts.et);
   t8->SetBranchAddress("eta",genParts.eta);
   t8->SetBranchAddress("phi",genParts.phi);
   t8->SetBranchAddress("mult",&genParts.n);

   cout<<"a"<<endl;

   for(int iev = 0; iev < t1->GetEntries(); ++iev){
      t1->GetEntry(iev);
      t2->GetEntry(iev);
      t3->GetEntry(iev);
      t4->GetEntry(iev);
      t5->GetEntry(iev);
      t6->GetEntry(iev);
      t7->GetEntry(iev);
      t8->GetEntry(iev);
      cout<<"a"<<endl;

      float EThbheRecHit = 0;

      for(int i = 0; i < hbheRecHit.n; ++i){
         float e = hbheRecHit.e[i];
         float et = hbheRecHit.et[i];
         float eta = hbheRecHit.eta[i];
         float phi = hbheRecHit.phi[i];
         EThbheRecHit += et;
      }

      float ETmyTowers = 0;
      int NmyTowers = 0;
      for(int i = 0; i < myTowers.n; ++i){
         float e = myTowers.e[i];
         float et = myTowers.et[i];
         float eta = myTowers.eta[i];
         float phi = myTowers.phi[i];

	 if(fabs(eta) < 2.87) continue;
	 ETmyTowers += et;

	 if(e < towerCut) continue;
	 NmyTowers += et;
      }

      float EThfRecHit = 0;
      int NhfRecHit = 0;

      for(int i = 0; i < hfRecHit.n; ++i){
         float e = hfRecHit.e[i];
         float et = hfRecHit.et[i];
         float eta = hfRecHit.eta[i];
         float phi = hfRecHit.phi[i];
         EThfRecHit += et;
         if(e < towerCut) continue;
         NhfRecHit++;
      }

      float ETgenParts = 0;
      int NgenParts = 0;

      for(int i = 0; i < genParts.n; ++i){
         float e = genParts.e[i];
         float et = genParts.et[i];
         float eta = genParts.eta[i];
         float phi = genParts.phi[i];

         if(fabs(eta) < 2.87) continue;
         NgenParts++;
	 ETgenParts += et;
      }

      nt->Fill(ETmyTowers,NmyTowers,EThfRecHit,NhfRecHit,NgenParts,ETgenParts);
   }


   cout<<"Main loop complete"<<endl;

   t1->AddFriend(nt);


   TCanvas* c1 = new TCanvas("c1","",1200,800);
   c1->Divide(5,3);

   c1->cd(1);
   t1->Draw(Form("NhfRecHit"),trig);
   c1->cd(2);
   //   t1->Draw("hf.eta:hf.phi",trig,"colz");
  c1->cd(3);
  //  t1->Draw("Sum$(hf.e):Npart",trig,"colz");
  c1->cd(4);
  t1->Draw("rho[1]+rho[9]:Npart",trig,"colz");
  c1->cd(5);
  t1->Draw("rho[6]:Npart",trig,"colz");
  c1->cd(6);
  t1->Draw("rho[1]+rho[9]+sigma[1]+sigma[9]:Npart",trig,"colz");

  c1->cd(7);
  t1->Draw("rho[6]+sigma[6]:Npart",trig,"colz");
  c1->cd(8);
  t1->Draw("hiEB:Npart",trig,"colz");

  c1->cd(9);
  t1->Draw("hiNpix:Npart",trig,"colz");

  c1->cd(10);
  t1->Draw(Form("NmyTowers>>ha0"));
  t1->Draw(Form("NmyTowers>>ha1"),trig,"same");

  c1->cd(11);
  t1->Draw("ETmyTowers:hiHF",trig,"colz");
  
  c1->cd(12);
  t1->Draw(Form("NmyTowers:NhfRecHit"),trig,"colz");


  if(MC){
  TCanvas * c2 = new TCanvas("c2","",900,600);
  c2->Divide(3,2);
  c2->cd(1);
  t1->Draw("hiHF:Npart>>h2o0",trig,"colz");

  c2->cd(2);

  t1->Draw("NgenParts:Npart>>h2a0",trig,"colz");
  c2->cd(3);
  t1->Draw("hiHF:NgenParts>>h2b0",trig,"colz");

  c2->cd(4);
  t1->Draw("ETgenParts:Npart>>h2c0",trig,"colz");
  c2->cd(5);
  t1->Draw("hiHF:ETgenParts>>h2d0",trig,"colz");

  c2->cd(6);
  t1->Draw("NgenParts>>hn0",trig);
  t1->Draw("ETgenParts>>he0",trig);

  }
  outf->Write();

  return 1;
}

