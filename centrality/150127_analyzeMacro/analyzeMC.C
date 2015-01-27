

void analyze(const char* infile = "test.root", const char* outfile = "outpu3.root"){

   const char* trig = "L1Tech_BSC_minBias_threshold1.v0";
   double towerCut = 3;

   bool MC = true;

   TFile * inf = new TFile(infile);
   TFile* outf = new TFile(outfile,"recreate");

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
   TTree* t3 = (TTree*)inf->Get("rechits/ee");
   TTree* t4 = (TTree*)inf->Get("rechits/bkg");
   TTree* t5 = (TTree*)inf->Get("rechits/tower");
   if(MC){
      TTree* t6 = (TTree*)inf->Get("mc/hi");
      t1->AddFriend(t6);
   }

   t1->AddFriend(t2);
   t1->AddFriend(t3);
   t1->AddFriend(t4);
   t1->AddFriend(t5);

   TCanvas* c1 = new TCanvas("c1","",1200,800);
   c1->Divide(5,3);

   c1->cd(1);
   t1->Draw(Form("Sum$(hf.e > %f)",towerCut),trig);
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
  t1->Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.e > %f)>>ha0",towerCut));
  t1->Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.e > %f)>>ha1",towerCut),trig,"same");

  c1->cd(11);
  t1->Draw("Sum$(tower.et * (abs(tower.eta) > 2.87)):hiHF",trig,"colz");
  
  c1->cd(12);
  t1->Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.e > %f):Sum$(hf.e > %f)",towerCut,towerCut),trig,"colz");


  if(MC){
  TCanvas * c2 = new TCanvas("c2","",900,600);
  c2->Divide(3,2);
  c2->cd(1);
  t1->Draw("hiHF:Npart>>h2o0",trig,"colz");

  c2->cd(2);

  t1->Draw("Sum$(hi.eta > 2.87 && hi.eta < 5.205):Npart>>h2a0",trig,"colz");
  c2->cd(3);
  t1->Draw("hiHF:Sum$(hi.eta > 2.87 && hi.eta < 5.205)>>h2b0",trig,"colz");

  c2->cd(4);
  t1->Draw("Sum$(hi.pt* (hi.eta > 2.87 && hi.eta < 5.205)):Npart>>h2c0",trig,"colz");
  c2->cd(5);
  t1->Draw("hiHF:Sum$(hi.pt* (hi.eta > 2.87 && hi.eta < 5.205))>>h2d0",trig,"colz");

  c2->cd(6);
  t1->Draw("Sum$(hi.eta > 2.87 && hi.eta < 5.205)>>hn0",trig);
  t1->Draw("Sum$(hi.pt* (hi.eta > 2.87 && hi.eta < 5.205))>>he0",trig);

  }
  outf->Write();

}

