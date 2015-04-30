///////////////////////////////////////////////////////////////////
// forest2yskim.C                                                //
// Creator : Yongsun Kim (MIT), jazzitup@mit.edu                 //
// Function : Transform hiForest files into yskim file           //
// yskims for MinBias1, Minbias2 and photon jet skims            //
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <TMath.h>
#include "../HiForestAnalysis/hiForest.h"
#include "../gammaJetAnalysis/CutAndBinCollection2012.h"
#include <time.h>
#include <TRandom3.h>

using namespace std;

static const int MAXTRK  = 10000;   // This is very enough.
static const int MAXMTRK = 30000;   // Again this is very enough for 10 mixing
static const int MAXMJET = 2000;
static const long MAXTREESIZE = 10000000000;


vector<jetKinem> nullVec;




void forest2yskim_jetSkim_forestV3(TString inputFile_="forestFiles/pA/pA_photonSkimForest_v85_skimPhotonPt50_eta1.5.root",
				   std::string MinbiasFname = "skim_trackJet_minbiasTrackJet_mc.root",
				   float cutphotonPt  = 35,  // default value dropped to 35GeV  for later photon energy smearing/scaling
				   std::string outname = "testPhotonSkim.root",
				   sampleType colli=kPADATA,
				   bool doMix = false,
				   bool doJetResCorrection = 1,  // = L2L3 * MC nonclosure correction  jet energy correction is done by default from Oct 19th (YS)
				   int smearingCentBin = -1, //0=0-10%, 1=10-30%, 2=30-50%, 3=50-100%, 4=0-30%, 5=30-100%  : Jet pT and phi smearing!
				   bool useGenJetColl = 0
  )
{
  bool isMC=true;
  if ((colli==kPPDATA)||(colli==kPADATA)||(colli==kHIDATA))
    isMC=false;

  int seconds = time(NULL);         cout << " time = " <<seconds%10000<< endl;
  TRandom3 rand(seconds%10000);
  TString datafname  = "";
  float cutphotonEta = 1.44;
  float preCutPhotonEt = 30;

  const int nMaxPho = 100;


  HiForest *c;
  if((colli==kPADATA)||(colli==kPAMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPPb, isMC );
  }
  else if  ((colli==kPPDATA)||(colli==kPPMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPP, isMC );
  }
  else if  ((colli==kHIDATA)||(colli==kHIMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPbPb, isMC );
    c->GetEnergyScaleTable("../photonEnergyScaleTable_lowPt_v6.root");
  }
  else {
    cout << " Error!  No such collision type" << endl;
    return;
  }

  c->InitTree();


  // vertex and centrality vtxCentWeighting
  TFile* fWeight = new TFile("../gammaJetAnalysis/vertexReweightingHistogram_pthatweighted.root");

  TH1D* hWeight_vtx_data_pp = (TH1D*)fWeight->Get("vertexHistoData_pp");
  TH1D* hWeight_vtx_mc_pp = (TH1D*)fWeight->Get("vertexHistoMC_pp");

  TH1D* hWeight_vtx_data_ppb = (TH1D*)fWeight->Get("vertexHistoData_ppb");
  TH1D* hWeight_vtx_mc_ppb = (TH1D*)fWeight->Get("vertexHistoMC_ppb");

  TH1D* hWeight_vtx_data_pbpb = (TH1D*)fWeight->Get("vertexHistoData_pbpb");
  TH1D* hWeight_vtx_mc_pbpb = (TH1D*)fWeight->Get("vertexHistoMC_pbpb");
  TH1D* hWeight_cent_data_pbpb = (TH1D*)fWeight->Get("centBinHistoData_pbpb");
  TH1D* hWeight_cent_mc_pbpb = (TH1D*)fWeight->Get("centBinHistoMC_pbpb");




  // L2L3 correction
  TFile* fL2L3pp = new TFile("../gammaJetAnalysis/corrL2L3/Casym_pp_double_hcalbins_algo_ak3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
  TH1D * c_etapp=(TH1D*)fL2L3pp->Get("C_asym");
  TF1* fptpp = new TF1("fptpp","1-[0]/pow(x,[1])",20,300);
  fptpp->SetParameters(0.06971,0.8167);

  TFile* fL2L3pA = new TFile("../gammaJetAnalysis/corrL2L3/Casym_pPb_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
  TH1D * c_etapA=(TH1D*)fL2L3pA->Get("C_asym");
  TF1* fptpA = new TF1("fptpA","1-[0]/pow(x,[1])",20,300);
  fptpA->SetParameters(0.3015, 0.8913);

  TFile* fL2L3Ap = new TFile("../gammaJetAnalysis/corrL2L3/Casym_Pbp_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
  TH1D * c_etaAp=(TH1D*)fL2L3Ap->Get("C_asym");
  TF1* fptAp = new TF1("fptAp","1-[0]/pow(x,[1])",20,300);
  fptAp->SetParameters(0.3015, 0.8913);

  // pA MC
  TF1 * fgaus=new TF1("fgaus_pA","gaus(0)",-20,20);
  fgaus->SetParameters(1,0,1);
  TF1 * fsmear_pA = new TF1("fsmear_pA","[0]/pow(x,[1])",50,300);
  fsmear_pA->SetParameters(1.052,0.5261);


  /*
     TCanvas* c11 = new TCanvas("c11","",1200,400);   // valiation of smearing factors
     c11->Divide(3,1);
     c11->cd(1);
     c_etapp->Draw();
     c11->cd(2);
     c_etaAp->Draw();
     c11->cd(3);
     c_etapA->Draw();
     c11->SaveAs("f1.gif");
  */


  // Create a new root file
  TFile* newfile_data = new TFile(outname.data(),"recreate");

  TTree* newtreePhoton;
  float newPt[nMaxPho];  // <<= temporary space
  int order[nMaxPho];
  float corrPt[nMaxPho];
  newtreePhoton = c->photonTree->CloneTree(0);
  newtreePhoton->SetName("yPhotonTree");
  newtreePhoton->SetMaxTreeSize(MAXTREESIZE);
  newtreePhoton->Branch("order",  order, "order[nPhotons]/I");
  newtreePhoton->Branch("corrPt", corrPt,"corrPt[nPhotons]/F");

  TTree* treeFullJet;
  if (   (colli==kPPDATA) || (colli==kPPMC) ) {
  //if (   (colli==kPPDATA) || (colli==kPPMC) ) {
    //treeFullJet = c->ak3jetTree->CloneTree(0);
    treeFullJet = c->ak3PFJetTree->CloneTree(0);
    cout << "pp collision.  Using ak3PF Jet Algo" << endl<<endl;
  }
  else {
    //treeFullJet = c->akPu3jetTree->CloneTree(0);
    treeFullJet = c->akPu3PFJetTree->CloneTree(0);
    cout << "PbPb or pPb collision. Using akPu3PF Jet Algo" << endl<<endl;
  }
  treeFullJet->SetName("fullJet");
  treeFullJet->SetMaxTreeSize(MAXTREESIZE);
#if 1
  TTree* treeGenp;
  if (   (colli==kHIMC ) || (colli==kPPMC) || (colli==kPAMC) ) {
    treeGenp =  c->genParticleTree->CloneTree(0);
    treeGenp->SetName("genparTree");
    //treeGenp =  c->genpTree->CloneTree(0);
    //treeGenp->SetName("genparTree");
    treeGenp->SetMaxTreeSize(MAXTREESIZE);
  }
#endif


  // jet tree!
  int nJet;
  const int MAXJET = 50000; // to accomodate 100 smeared jets, need to be careful with ram
  float jetPt[MAXJET];
  float jetEta[MAXJET];
  float jetPhi[MAXJET];
  float jetDphi[MAXJET];
  int jetSubid[MAXJET];
  float jetRefPt[MAXJET];
  float jetRefEta[MAXJET];
  float jetRefPhi[MAXJET];
  float jetRefDphi[MAXJET];
  float jetRefPartonPt[MAXJET];
  int  jetRefPartonFlv[MAXJET];

  TTree *newtreeJet = new TTree("yJet","jets");
  newtreeJet->SetMaxTreeSize(MAXTREESIZE);
  newtreeJet->Branch("nJet",&nJet,"nJet/I");
  newtreeJet->Branch("pt",jetPt,"pt[nJet]/F");
  newtreeJet->Branch("eta",jetEta,"eta[nJet]/F");
  newtreeJet->Branch("phi",jetPhi,"phi[nJet]/F");
  newtreeJet->Branch("dphi",jetDphi,"dphi[nJet]/F");
  if ( isMC )  {
    newtreeJet->Branch("subid",jetSubid,"subid[nJet]/I");
    newtreeJet->Branch("refPt",jetRefPt,"refPt[nJet]/F");
    newtreeJet->Branch("refEta",jetRefEta,"refEta[nJet]/F");
    newtreeJet->Branch("refPhi",jetRefPhi,"refPhi[nJet]/F");
    newtreeJet->Branch("refDphi",jetRefDphi,"refDphi[nJet]/F");
    newtreeJet->Branch("refPartonPt",jetRefPartonPt,"refPartonPt[nJet]/F");
    newtreeJet->Branch("refPartonFlv",jetRefPartonFlv,"refPartonFlv[nJet]/I");
  }



  int nMjet;
  float mJetPt[MAXMJET];
  float mJetEta[MAXMJET];
  float mJetPhi[MAXMJET];
  float mJetDphi[MAXMJET];

  TTree * tmixJet = new TTree("mJet","Jet from minbias events");
  tmixJet->SetMaxTreeSize(MAXTREESIZE);
  tmixJet->Branch("nJet",&nMjet,"nJet/I");
  tmixJet->Branch("pt",mJetPt,"pt[nJet]/F");
  tmixJet->Branch("eta",mJetEta,"eta[nJet]/F");
  tmixJet->Branch("phi",mJetPhi,"phi[nJet]/F");
  tmixJet->Branch("dphi", mJetDphi, "dphi[nJet]/F");



  // Imb = Input MinBias events
  EvtSel          evtImb;
  Int_t           nJetImb;
  Float_t         jetPtImb[100];
  Float_t         jetEtaImb[100];
  Float_t         jetPhiImb[100];
  TBranch        *b_evt;
  TBranch        *b_nJetImb;
  TBranch        *b_jetPtImb;
  TBranch        *b_jetEtaImb;
  TBranch        *b_jetPhiImb;


  int nCentBins =  nCentBinSkim;
  if ((colli==kPADATA)||(colli==kPAMC)) {
    nCentBins = nCentBinSkimPA;
  }

  TChain   *tjmb[100][nVtxBin+1];
  int nMB[100][nVtxBin+1] ; //= 199109;
  int mbItr[100][nVtxBin+1];
  if ( doMix ) {
    cout <<"  Tree initialization for MinBias mixing" << endl;
    for( int icent = 0 ; icent< nCentBins ; icent++) {
      for( int ivz = 1 ; ivz<=nVtxBin ; ivz++) {

	tjmb[icent][ivz] = new TChain(Form("trkAndJets_first_cBin2icent%d_ivz%d",icent,ivz));
	tjmb[icent][ivz]->Add(MinbiasFname.data());
	tjmb[icent][ivz]->SetBranchAddress("evt", &evtImb,&b_evt);
	tjmb[icent][ivz]->SetBranchAddress("nJet",   &nJetImb,   &b_nJetImb);
	tjmb[icent][ivz]->SetBranchAddress("jetPt",  &jetPtImb,  &b_jetPtImb);
	tjmb[icent][ivz]->SetBranchAddress("jetEta", &jetEtaImb, &b_jetEtaImb);
	tjmb[icent][ivz]->SetBranchAddress("jetPhi", &jetPhiImb, &b_jetPhiImb);


	nMB[icent][ivz] = tjmb[icent][ivz]->GetEntries();
	cout << "number of evetns in (icent = " << icent << ", ivtxZ = "<< ivz << ")  = " << nMB[icent][ivz] << endl;
	int primeSeed = rand.Integer(37324);
	mbItr[icent][ivz] = primeSeed%(nMB[icent][ivz]);
	cout <<" initial itr = " << mbItr[icent][ivz] << endl;
      }
    }
  }
  else
    cout << endl << endl << "  Mixing process is skipped" << endl << endl << endl ;


  // reweighting factor should go here
  int eTot(0), eSel(0);
  EvtSel evt;
  GammaJet gj;
  Isolation isol;
  TTree *tgj;
  tgj = new TTree("tgj","gamma jet tree");
  tgj->SetMaxTreeSize(MAXTREESIZE);
  tgj->Branch("evt",&evt.run,"run/I:evt:cBin:pBin:trig/O:offlSel:noiseFilt:anaEvtSel:vz/F:vtxCentWeight/F:hf4Pos:hf4Neg:hf4Sum:ptHat:ptHatWeight");
  tgj->Branch("lpho",&gj.photonEt,"photonEt/F:photonRawEt:photonEta:photonPhi:hovere:r9:sigmaIetaIeta:sumIso:genIso:genPhotonEt:genMomId/I:lJetPt/F:lJetEta:lJetPhi:lJetDphi:lJetSubid/I");
  tgj->Branch("isolation",&isol.cc1,"cc1:cc2:cc3:cc4:cc5:cr1:cr2:cr3:cr4:cr5:ct1PtCut20:ct2PtCut20:ct3PtCut20:ct4PtCut20:ct5PtCut20:ecalIso:hcalIso:trackIso");  // ecalIso,hcalIso,trackIso are the pp style isolation


  float vzCut = vtxCutPhotonAna;
  TH1F* hvz = new TH1F("hvz","",nVtxBin,-vzCut,vzCut);
  // event plane hitogram
  TH1F* hEvtPlnBin = new TH1F("hEvtPlnBin", "", nPlnBin, -PI/2., PI/2.);
  // jet algos
  Jets* theJet;
  Jets* genJetTree;
  if (   (colli==kPPDATA) || (colli==kPPMC) ) {
    theJet = &(c->ak3PF) ;
    cout << "pp collision.  Using ak3PF Jet Algo" << endl<<endl;
  }
  else {
    theJet = &(c->akPu3PF) ;
    cout << "PbPb or pPb collision. Using akPu3PF Jet Algo" << endl<<endl;
  }
  genJetTree = &(c->akPu3PF);


  // Loop starts.
  int nentries = c->GetEntries();
  cout << "number of entries = " << nentries << endl;
  for (Long64_t jentry = 0 ; jentry < nentries; jentry++) {
    eTot++;
    if (jentry% 2000 == 0)  {
      cout <<jentry<<" / "<<nentries<<" "<<setprecision(2)<<(double)jentry/nentries*100<<endl;
    }

    c->GetEntry(jentry);
    // Select events with a generated photon in mid-rapidity
    bool genPhotonFlag=false;
    if ( !isMC )   // fixed the most stupid error
      genPhotonFlag = true;
    else {
#if 0 // there is no genp tree, but is higentree
      for ( int g=0 ; g< c->genp.nPar ; g++) {
	if ( c->genp.id[g] != 22 )
	  continue;
	if ( fabs( c->genp.momId[g] ) > 22 )
	  continue;
	if ( fabs( c->genp.status[g] ) != 1 )
	  continue;
	if ( fabs( c->genp.eta[g] ) > cutphotonEta )
	  continue;
	if ( c->genp.et[g] < 35 )
	  continue;
#endif
     for ( int g=0 ; g< c->genparticle.mult ; g++) {
	if ( c->genparticle.pdg[g] != 22 )
	  continue;
//	if ( fabs( c->genparticle.momId[g] ) > 22 )
//	  continue;
//	if ( fabs( c->genparticle.status[g] ) != 1 )
//	  continue;
	if ( fabs( c->genparticle.eta[g] ) > cutphotonEta )
	  continue;
	if ( c->genparticle.pt[g] < cutphotonPt )
	  continue;
	genPhotonFlag = true;
      }
    }
    if ( !genPhotonFlag)
      continue;

    evt.clear();
    evt.run   = c->evt.run;
    evt.evt = c->evt.evt;
    evt.hf4Pos = c->evt.hiHFplusEta4;
    evt.hf4Neg = c->evt.hiHFminusEta4;
    evt.hf4Sum = evt.hf4Pos + evt.hf4Neg;
    evt.cBin = -99;
    evt.pBin   = -99 ;
    if ((colli==kHIDATA)||(colli==kHIMC))   {
      evt.cBin = c->evt.hiBin;
      evt.pBin   = hEvtPlnBin->FindBin( c->evt.hiEvtPlanes[theEvtPlNumber] ) ;
    }
    else if ((colli==kPADATA)||(colli==kPAMC))   {
      evt.cBin =  getHfBin(evt.hf4Sum);
      // if (  ((evt.cBin) < 0) || (evt.cBin) > 18 )
      //   cout << " Check the pA centrality..  cbin = " << evt.cBin << endl;
    }

    evt.trig = 0;
    evt.offlSel = (c->skim.pcollisionEventSelection > 0);
    evt.noiseFilt = (c->skim.pHBHENoiseFilter > 0);
    evt.anaEvtSel = c->selectEvent() && evt.trig;
    evt.vz = c->evt.vz;

    int cBin = evt.cBin;
    int vzBin = hvz->FindBin(evt.vz)  ;
    hvz->Fill(evt.vz) ;

	// this was the problem!!! all c->selectEvent() are 0
    //if ( ( (colli==kHIDATA)||(colli==kHIMC)||(colli==kPADATA)||(colli==kPAMC) || (colli==kPPMC) ) && ( c->selectEvent() == 0 ))
    if ( ( (colli==kHIDATA)||(colli==kHIMC) ) && ( c->skim.pcollisionEventSelection == 0 ))
	continue;
    if ( ( (colli==kPADATA)||(colli==kPAMC) ) && ( c->skim.pPAcollisionEventSelectionPA == 0 ))
  	continue;
    if ( ( (colli==kPPMC) ) && ( c->skim.pcollisionEventSelection == 0 ))
      continue;
    if ( ( (colli==kPADATA)||(colli==kPPDATA) ) && ( c->skim.pVertexFilterCutGplus ==0 ) ) // No Pile up events
      continue;
    if ( (vzBin<1) || ( vzBin > nVtxBin) )
      continue;

    eSel++;      // OK.  This event is a collisional and no-noise event.


    // Reweight for vertex and centrality of MC
    evt.vtxCentWeight = 1;
    double wVtx=1;
    double wCent=1;
    if (colli ==kHIMC) {
      int vBin =  hWeight_vtx_data_pbpb->FindBin(evt.vz);
      wVtx =   hWeight_vtx_data_pbpb->GetBinContent(vBin) / hWeight_vtx_mc_pbpb->GetBinContent(vBin) ;
      wCent =  hWeight_cent_data_pbpb->GetBinContent(evt.cBin+1) / hWeight_cent_mc_pbpb->GetBinContent(evt.cBin+1) ;
    }
    else if ( colli ==kPPMC) {
      int vBin =  hWeight_vtx_data_pp->FindBin(evt.vz);
      wVtx =   hWeight_vtx_data_pp->GetBinContent(vBin) / hWeight_vtx_mc_pp->GetBinContent(vBin) ;
    }
    else if ( colli ==kPAMC) {
      int vBin =  hWeight_vtx_data_ppb->FindBin(evt.vz);
      wVtx =   hWeight_vtx_data_ppb->GetBinContent(vBin) / hWeight_vtx_mc_ppb->GetBinContent(vBin) ;
    }
    evt.vtxCentWeight = wVtx * wCent;

    evt.ptHat = -1;
    evt.ptHatWeight = 1;
    evt.ptHat = c->photon.ptHat;

   // if( colli == kPAMC && evt.ptHat > maxpthat ) // pA samples don't use ptHatCutter.C.
   //	continue;
    if (colli ==kHIMC) {
      if ( evt.ptHat < 50  )       evt.ptHatWeight = 9008/16237. ;
      else if ( evt.ptHat < 80  )       evt.ptHatWeight = 3750/85438. ;
      else   evt.ptHatWeight =  1191/140432. ;
    }
    else if ( colli == kPPMC) { // pp has only 4 pthat samples
      if ( evt.ptHat < 50  )       evt.ptHatWeight = 9008/9008. ;
      else if ( evt.ptHat < 80  )       evt.ptHatWeight = 3750/40109. ;
      else   evt.ptHatWeight = 1191/66934. ;
    }
    else if ( colli == kPAMC) { 
	if ( evt.ptHat > 30 && evt.ptHat < 50  )       evt.ptHatWeight = 62744/62744. ;
        else if ( evt.ptHat > 50 && evt.ptHat < 80  )       evt.ptHatWeight = 29499/107309. ;
        else if ( evt.ptHat > 80 && evt.ptHat < 120  )       evt.ptHatWeight = 7640/106817. ;
        else if ( evt.ptHat > 120 && evt.ptHat < 170  )       evt.ptHatWeight = 1868/104443. ;
        else   evt.ptHatWeight = 649/139647. ;
    }



   for (int j=0;j< c->photon.nPhotons;j++) {

      if (  ( c->photon.pt[j] > preCutPhotonEt ) && ( fabs( c->photon.eta[j] ) < cutphotonEta ) ) {
	newPt[j] = c->getCorrEt(j);
      }
      else
	newPt[j] = c->photon.pt[j] - 10000;

      if ( (c->isSpike(j)) || (c->photon.hadronicOverEm[j]>0.2) ) //||  (c->photon.isEle[j]) )
	newPt[j] = newPt[j] - 20000;
      if (c->photon.seedTime[j] ==0 )   // clustering bug
	newPt[j] = newPt[j] - 30000;

      corrPt[j] = newPt[j];
    }

    TMath::Sort(c->photon.nPhotons, newPt, order);


    // Select the leading photon
    gj.clear();
    int leadingIndex=-1;
    for (int j=0;j<c->photon.nPhotons;j++) {
      if ( c->photon.pt[j]  < preCutPhotonEt ) continue;
      if ( fabs(c->photon.eta[j]) > cutphotonEta ) continue;
      if (c->isSpike(j)) continue;
      //   if (!(c->isLoosePhoton(j))) continue;
      if (c->photon.hadronicOverEm[j]>0.1) continue;
      if ((c->photon.rawEnergy[j]/c->photon.energy[j])<0.5) continue;



      // sort using corrected photon pt
      float theCorrPt= corrPt[j];
      if ( theCorrPt > gj.photonEt) {
	gj.photonEt = theCorrPt;
	leadingIndex = j;
      }
    }
    //    if ( (gj.photonEt < cutphotonPt) )    <== This cut ruins the ptHat weighting factor
    //  continue;


    /// Save leading photons
    if (leadingIndex!=-1) {
      gj.photonRawEt=c->photon.pt[leadingIndex];
      gj.photonEta=c->photon.eta[leadingIndex];
      gj.photonPhi=c->photon.phi[leadingIndex];
      gj.hovere=c->photon.hadronicOverEm[leadingIndex];
      gj.r9=c->photon.r9[leadingIndex];
      gj.sigmaIetaIeta=c->photon.sigmaIetaIeta[leadingIndex];
      gj.sumIsol = (c->photon.cr4[leadingIndex]+c->photon.cc4[leadingIndex]+c->photon.ct4PtCut20[leadingIndex]) / 0.9;
      gj.genIso   = c->photon.genCalIsoDR04[leadingIndex];
      gj.genPhotonEt = c->photon.genMatchedPt[leadingIndex];
      gj.genMomId = c->photon.genMomId[leadingIndex];
      isol.Set(c,leadingIndex);
    }
    else {
      gj.clear();
    }
    if ( (colli==kPADATA) && ( evt.run > 211256 ) )  {
      gj.photonEta = - gj.photonEta;
    }



    ///////////////////// Skimmed Jet tree ///////////////////////////////////
    nJet = 0 ;

    int jetEntries = 0;
    if (useGenJetColl )    jetEntries = theJet->ngen;
    else                   jetEntries = theJet->nref;

    int nSmear = 1;
    if(smearingCentBin != -1)
      nSmear = 100;

      for(int iSmear =0; iSmear < nSmear; iSmear++){ // iSmear loop have to be here, before Jet loop. mis-ordered for loops ruins jetrefpt values.
    for (int ij=0; ij< jetEntries ; ij++) {
      if ( gj.photonEt < 0 )    continue ;    // If there is no photon in this event
	if (  useGenJetColl )   {
	  jetPt[nJet] = theJet->genpt[ij];
	  jetEta[nJet] = theJet->geneta[ij];
	  jetPhi[nJet] = theJet->genphi[ij];
	}
	else  {
	  jetPt[nJet] = theJet->jtpt[ij];
	  jetEta[nJet] = theJet->jteta[ij];
	  jetPhi[nJet] = theJet->jtphi[ij];
	}

	// Smear phi
	Double_t newPhi = jetPhi[nJet] ;
	if( smearingCentBin != -1 )
	{
	  Double_t phiSmear  = TMath::Sqrt((cphi_pbpb[smearingCentBin]*cphi_pbpb[smearingCentBin] - cphi_pp*cphi_pp)
					   + (sphi_pbpb[smearingCentBin]*sphi_pbpb[smearingCentBin] - sphi_pp*sphi_pp)/jetPt[nJet]
					   + (nphi_pbpb[smearingCentBin]*nphi_pbpb[smearingCentBin] - nphi_pp*nphi_pp)/(jetPt[nJet]*jetPt[nJet]));
	  newPhi  =  jetPhi[nJet] +   rand.Gaus(0, phiSmear);
	  while ( fabs(newPhi) > PI )  {
	    if ( newPhi > PI )  newPhi = newPhi - 2*PI;
	    if ( newPhi < -PI )  newPhi = newPhi + 2*PI;
	  }
	}
	jetPhi[nJet] = newPhi;


	// smear the jet pT
	//float smeared = jetPt[nJet] * rand.Gaus(1,addJetEnergyRes/jetPt[nJet])   *  rand.Gaus(1, addFlatJetEnergyRes) ;
	Double_t smeared = jetPt[nJet];
	if( smearingCentBin != -1 )
	{
	  Double_t smearSigma = TMath::Sqrt((c_pbpb[smearingCentBin]*c_pbpb[smearingCentBin] - c_pp*c_pp)
					    + (s_pbpb[smearingCentBin]*s_pbpb[smearingCentBin] - s_pp*s_pp)/jetPt[nJet]
					    + (n_pbpb[smearingCentBin]*n_pbpb[smearingCentBin] - n_pp*n_pp)/(jetPt[nJet]*jetPt[nJet]));
	  smeared = jetPt[nJet] * rand.Gaus(1, smearSigma);
	}
	// then multiply jet energy sclae

	// resCorrection
	float resCorrection =1. ;
	float l2l3Corr =1 ;

	if  (doJetResCorrection)   {
	  // L2L3 correction!
	  if ( colli == kPPDATA)   {
	    l2l3Corr = c_etapp->GetBinContent(c_etapp->FindBin(jetEta[nJet])) * fptpp->Eval( jetPt[nJet]);
	  }
	  else if ( colli == kPADATA)   {
	    if ( evt.run > 211256 )
	      l2l3Corr = c_etapA->GetBinContent(c_etapA->FindBin(jetEta[nJet]))  * fptpA->Eval( jetPt[nJet]);
	    else
	      l2l3Corr = c_etaAp->GetBinContent(c_etaAp->FindBin(jetEta[nJet]))  * fptAp->Eval( jetPt[nJet]);
	  }
	  else if ( colli == kPAMC)
	    l2l3Corr = 1 + (fsmear_pA->Eval( jetPt[nJet] )) * fgaus->GetRandom()  ;

	  // do the residual correction
	  if ((colli==kHIDATA)||(colli==kHIMC))  { 
	    if ( evt.cBin  < 12 )   // central
	      resCorrection  =  1.04503 -1.6122  /(sqrt(jetPt[nJet])) + 9.27212 / (jetPt[nJet]);  //1.04503    -1.6122    9.27212
	    else                  // peripheral
	      resCorrection  =  1.00596 -0.653191/(sqrt(jetPt[nJet])) + 4.35373 / (jetPt[nJet]);  //1.00596     -0.653191  4.35373
	  }
	  else if ((colli==kPPDATA)||(colli==kPPMC)){  // do the residual correction
	    resCorrection  = 0.993609  +0.158418/(sqrt(jetPt[nJet])) + 0.335479 / (jetPt[nJet]);//	  0.993609   0.158418   0.335479
	  }
	  else if ((colli==kPADATA)||(colli==kPAMC)){
	    resCorrection = 0.997738 + 0.0221806/(sqrt(jetPt[nJet])) - 0.877999/ (jetPt[nJet]); //C : 0.997738, S : 0.0221806, N : -0.877999 //function derived as a function of gen pt 
	    //resCorrection = 0.981365 + 0.342746/(sqrt(jetPt[nJet])) - 2.63018 / (jetPt[nJet]); //C : 0.981365, S : 0.342746, N : -2.63018 //function derived as a function of gen pt //derived from l2l3 corrected sample. 
	    //resCorrection = 0.745753 + 6.91646/(sqrt(jetPt[nJet])) - 33.0167 / (jetPt[nJet]); //C : 0.745753, S : 6.91646, N : -33.0167 //function derived as a function of reco pt 
	  }
	
	} // doJetResCorrection

	// reflect eta!
	if ( (colli==kPADATA) && ( evt.run > 211256 ) )  {
	  jetEta[nJet] = -jetEta[nJet];
	}

	jetPt[nJet] = smeared * l2l3Corr /resCorrection;



	if ( jetPt[nJet] < cutjetPtSkim)  // double cutjetPtSkim = 15; Oct 19th
	  continue;
	if ( fabs( jetEta[nJet] ) > cutjetEtaSkim )     // double cutjetEtaSkim = 3.0; Oct 19th
	  continue;
	if ( getDR( jetEta[nJet], jetPhi[nJet], gj.photonEta, gj.photonPhi) < 0.5 )
	  continue;


	if (jetPt[nJet] >0)
	  jetDphi[nJet] = getAbsDphi( jetPhi[nJet], gj.photonPhi) ;
	else
	  jetDphi[nJet] = -1;


	if (  useGenJetColl )   {
	  jetSubid[nJet] = -9999;
	  jetRefPt[nJet] = -9999;
	  jetRefEta[nJet] = -9999;
	  jetRefPhi[nJet] = -9999;
	  jetRefPt[nJet] =  -9999;
	  jetRefPartonPt[nJet] = -9999;
	  jetRefPartonFlv[nJet] = -9999;
	}
	else {
	  jetRefPt[nJet] = theJet->refpt[ij];
	  jetRefEta[nJet] = theJet->refeta[ij];
	  jetRefPhi[nJet] = theJet->refphi[ij];
	  if (jetRefPt[nJet] >0)
	    jetRefDphi[nJet] = getAbsDphi( jetRefPhi[nJet] , gj.photonPhi) ;
	  else
	    jetRefDphi[nJet] = -1;
	  jetRefPartonPt[nJet] = theJet->refparton_pt[ij];
	  jetRefPartonFlv[nJet] = theJet->refparton_flavor[ij];

	  jetSubid[nJet] = -9999;
	  if (jetRefPt[nJet] >0)  {   // Find the collisional subid of this gen jet!!
	    for ( int igen=0; igen < genJetTree->ngen ; igen++) {
	      if ( jetRefPt[nJet] == genJetTree->genpt[igen] )
		jetSubid[nJet] = genJetTree->gensubid[igen] ;
	    }

	    if ( jetSubid[nJet] == -9999 ) // No genJet matched!
	      cout << " WARNING!  This reco jet was not matched to anyone in the gen jet collection!!! " << endl;

	  }
	}

      nJet++ ;
      }// ij for loop
    }// iSmear for loop


    //////// Leading jet kinematics in dphi>7pi/8
    float maxJpt = 0;
    int jetLeadingIndex = -1;

    for (int ij=0; ij< nJet ; ij++) {
      if ( jetDphi[ij] < awayRange )  // const float awayRange= PI * 7./8.;
	continue;
      if ( fabs( jetEta[ij] ) > cutjetEta )  // double cutjetEta = 1.6;
	continue;
      if ( jetPt[ij] > maxJpt) {
	maxJpt = jetPt[ij] ;
	jetLeadingIndex = ij;
      }
    }
    if ( jetLeadingIndex > -1 ) {
      gj.lJetPt = jetPt[jetLeadingIndex];
      gj.lJetEta = jetEta[jetLeadingIndex];
      gj.lJetPhi = jetPhi[jetLeadingIndex];
      gj.lJetDphi =  jetDphi[jetLeadingIndex];
      gj.lJetSubid=  jetSubid[jetLeadingIndex];
    }
    else {
      gj.lJetPt = -1;
      gj.lJetEta = 999;
      gj.lJetPhi = 999;
      gj.lJetDphi = 0;
      gj.lJetSubid=  -99;
    }


    int nMixing = nMixing1;
    nMjet = 0;
    bool noSuchEvent = false;
    int iMix=0;
    int loopCounter=0;

    if ( !doMix )
      iMix = nMixing+1;   // Mixing step will be skipped

    while (iMix<nMixing)  {
      loopCounter++;
      if ( loopCounter > nMB[cBin][vzBin]+1) {
	iMix = 999999 ;
	noSuchEvent = true;
	cout << " no such event!! :  icent = " << cBin << ",  vzBin = " << vzBin << ",  pBin = " << evt.pBin << endl;
	continue;
      }

      mbItr[cBin][vzBin] = mbItr[cBin][vzBin] + 1;
      if ( mbItr[cBin][vzBin] == nMB[cBin][vzBin] )
	mbItr[cBin][vzBin] =  mbItr[cBin][vzBin] - nMB[cBin][vzBin];

      /// Load the minbias tracks!!
      tjmb[cBin][vzBin]->GetEntry(mbItr[cBin][vzBin]);



      // ok found the event!! ///////////
      loopCounter =0;  // Re-initiate loopCounter
      // Jet mixing
      for (int it = 0 ; it < nJetImb ; it++) {
	if ( gj.photonEt < 0 )  continue;
	// Smear phi
	Double_t newPhi = jetPhiImb[it];
	if( smearingCentBin != -1 )
	{
	  Double_t phiSmear  = TMath::Sqrt((cphi_pbpb[smearingCentBin]*cphi_pbpb[smearingCentBin] - cphi_pp*cphi_pp)
					   + (sphi_pbpb[smearingCentBin]*sphi_pbpb[smearingCentBin] - sphi_pp*sphi_pp)/jetPtImb[it]
					   + (nphi_pbpb[smearingCentBin]*nphi_pbpb[smearingCentBin] - nphi_pp*nphi_pp)/(jetPtImb[it]*jetPtImb[it]));
	  newPhi  =  jetPhiImb[it] +   rand.Gaus(0, phiSmear);
	  while ( fabs(newPhi) > PI )  {
	    if ( newPhi > PI )  newPhi = newPhi - 2*PI;
	    if ( newPhi < -PI )  newPhi = newPhi + 2*PI;
	  }
	}
	jetPhiImb[it] = newPhi;


	// smear the jet pT
	//float smeared = jetPtImb[it] * rand.Gaus(1,addJetEnergyRes/jetPtImb[it]) *  rand.Gaus(1, addFlatJetEnergyRes) ;
	Double_t smeared = jetPtImb[it];
	if( smearingCentBin != -1 )
	{
	  Double_t smearSigma = TMath::Sqrt((c_pbpb[smearingCentBin]*c_pbpb[smearingCentBin] - c_pp*c_pp)
					    + (s_pbpb[smearingCentBin]*s_pbpb[smearingCentBin] - s_pp*s_pp)/jetPtImb[it]
					    + (n_pbpb[smearingCentBin]*n_pbpb[smearingCentBin] - n_pp*n_pp)/(jetPtImb[it]*jetPtImb[it]));
	  smeared = jetPtImb[it] * rand.Gaus(1, smearSigma);
	}
	float resCorrection =1. ;
	float l2l3Corr =1 ;

	if  (doJetResCorrection)   {
	  // L2L3
	  if ( colli == kPPDATA)   {
	    l2l3Corr = c_etapp->GetBinContent(c_etapp->FindBin(jetEtaImb[it])) * fptpp->Eval( jetPtImb[it]);
	  }
	  else if ( colli == kPADATA)   {
	    if ( evt.run > 211256 )
	      l2l3Corr = c_etapA->GetBinContent(c_etapA->FindBin(jetEtaImb[it]))  * fptpA->Eval( jetPtImb[it]);
	    else
	      l2l3Corr = c_etaAp->GetBinContent(c_etaAp->FindBin(jetEtaImb[it]))  * fptAp->Eval( jetPtImb[it]);
	  }
	  else if ( colli == kPAMC)
	    l2l3Corr = 1 + (fsmear_pA->Eval( jetPtImb[it] )) * fgaus->GetRandom()  ;
	
	  // Correction from MC closure
	  if ((colli==kHIDATA)||(colli==kHIMC))  { // do the residual correction
	    if ( evt.cBin  < 12 )   // central
	      resCorrection  =  1.04503 -1.6122  /(sqrt(jetPtImb[it])) + 9.27212 / (jetPtImb[it]);  //1.04503    -1.6122    9.27212
	    else                  // peripheral
	      resCorrection  =  1.00596 -0.653191/(sqrt(jetPtImb[it])) + 4.35373 / (jetPtImb[it]);  //1.00596     -0.653191  4.35373
	  }
	  else if ((colli==kPPDATA)||(colli==kPPMC)){  // do the residual correction
	    resCorrection  = 0.993609  +0.158418/(sqrt(jetPtImb[it])) + 0.335479 / (jetPtImb[it]);//          0.993609   0.158418   0.335479
	  }
	  else if ((colli==kPADATA)||(colli==kPAMC)){
	    resCorrection = 0.997738 + 0.0221806/(sqrt(jetPtImb[it])) - 0.877999 / (jetPtImb[it]); //C : 0.997738, S : 0.0221806, N : -0.877999 //function derived as a function of gen pt
	    //resCorrection = 0.981365 + 0.342746/(sqrt(jetPtImb[it])) - 2.63018 / (jetPtImb[it]); //C : 0.981365, S : 0.342746, N : -2.63018 //function derived as a function of gen pt//derived from l2l3 corrected sample.
	    //resCorrection = 0.745753 + 6.91646/(sqrt(jetPtImb[it])) - 33.0167 / (jetPtImb[it]); //C : 0.745753, S : 6.91646, N : -33.0167//function derived as a function of reco pt
	  }
	
	}


	float smearedCorrected  = smeared *l2l3Corr / resCorrection; // residual correction

	if ( smearedCorrected < cutjetPtSkim )  // double cutjetPtSkim = 15; Oct 19th
	  continue;
	if ( fabs( jetEtaImb[it] ) > cutjetEtaSkim )   // double cutjetEtaSkim = 3.0; Oct 19th
	  continue;
	if ( getDR( jetEtaImb[it], jetPhiImb[it], gj.photonEta, gj.photonPhi) < 0.5 )  // This cut added for consistency ; Oct 19th
	  continue;

	mJetPt[nMjet]    = smearedCorrected;
	mJetEta[nMjet]   = jetEtaImb[it];
	mJetPhi[nMjet]   = jetPhiImb[it];
	if  ( mJetPt[nMjet]>0 )
	  mJetDphi[nMjet]  =   getAbsDphi(mJetPhi[nMjet], gj.photonPhi) ;
	else
	  mJetDphi[nMjet]=-1;

	nMjet++; // < == Important!
      }


      iMix++;
    }
    if ( noSuchEvent )
      continue;

    tgj->Fill();
    newtreeJet->Fill();
    tmixJet->Fill();
    newtreePhoton->Fill();
    treeFullJet->Fill();
    if (   (colli==kHIMC ) || (colli==kPPMC) || (colli==kPAMC) )
      treeGenp->Fill();
  }



  newfile_data->Write();
  //   newfile_data->Close();   // <<=== If there is close() function. writing stucks in the middle of looping.. I don't know why!!
  cout << " Done! "<< endl;
  cout << "    " << eSel<<" out of total "<<eTot<<" events were analyzed."<<endl;
}
