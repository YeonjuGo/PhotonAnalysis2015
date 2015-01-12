// Author : Yeonju Go
// cut gen photon threshold from 25 to 35 GeV
// not including leading photon awayside jet condition

#include "TROOT.h"
#include "TString.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TDirectory.h"
#include <stdlib.h>
using namespace std;

//! pp 2.76 TeV PYTHIA cross sections
const double xsec[12][3] ={{2.034e-01,15,30}, //! 15
	{1.075e-02,30,50}, //! 30
	{1.025e-03,50,80}, //! 50
	{9.865e-05,80,120}, //! 80
	{1.129e-05,120,170}, //! 120
	{1.465e-06,170,220}, //! 170
	{2.837e-07,220,280}, //! 220
	{5.323e-08,280,370}, //! 280
	{5.934e-09,370,460}, //! 370
	{8.125e-10,460,540}, //! 460
	{1.467e-10,540,9999}, //! 540
	{0,9999,9999}
};

float GetXsec(float/*maxpthat*/);

void HI_weight_pPb_changeThreshold(const char *infile="/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT15_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root",
		const char *outfile="test_corr_MC.root",
		float maxpthat=30.,
		float xSection=1.075e-02
		)
{
	TFile *fin = TFile::Open(infile, "readonly");
	cout<<infile<<endl;

	cout<<endl;
	cout<<endl;
	fin->ls();
	cout<<endl;
	cout<<endl;

	const int ndir=6;

	string DirName[ndir]= {
	//	"ak3CaloJetAnalyzer",
		"ak3PFJetAnalyzer",
	//	"ak4CaloJetAnalyzer",
		"ak4PFJetAnalyzer",
	//	"ak5CaloJetAnalyzer",
		"ak5PFJetAnalyzer",
	//	"akPu3CaloJetAnalyzer",
		"akPu3PFJetAnalyzer",
	//	"akPu4CaloJetAnalyzer",
		"akPu4PFJetAnalyzer",
	//	"akPu5CaloJetAnalyzer",
		"akPu5PFJetAnalyzer",
	};

	TFile *fout=new TFile(outfile,"RECREATE");

	TTree *tr_in=0, *tr_out=0, *tr_pho=0;
	Int_t           hiBin;
	Float_t           hiHF;

	TTree *tr_ev = 0;
	tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
	tr_ev->SetBranchAddress("hiHF",&hiHF);
	tr_ev->SetBranchAddress("hiBin",&hiBin);
	tr_ev->SetBranchStatus("*",0,0);
	tr_ev->SetBranchStatus("hiBin",1,0);
	tr_ev->SetBranchStatus("hiHF",1,0);

	Int_t	pPAcollisionEventSelectionPA, pPAprimaryVertexFilter;
	TTree *tr_skim = 0;
	tr_skim = (TTree*)fin->Get("skimanalysis/HltTree");
	tr_skim->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
	tr_skim->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
	tr_skim->SetBranchStatus("*",0,0);
	tr_skim->SetBranchStatus("pPAcollisionEventSelectionPA",1,0);
	tr_skim->SetBranchStatus("pPAprimaryVertexFilter",1,0);


        float ptgen[1000]; 
	float etagen[1000];
	float phigen[1000];
        int mult;
	int pdg[1000];
	TTree *tr_gen = 0;
        tr_gen = (TTree*)fin->Get("HiGenParticleAna/hi");
        tr_gen->SetBranchAddress("mult", &mult);
        tr_gen->SetBranchAddress("pdg", pdg);
        tr_gen->SetBranchAddress("pt", ptgen);
        tr_gen->SetBranchAddress("eta", etagen);
	tr_gen->SetBranchStatus("*",0,0);
	tr_gen->SetBranchStatus("mult",1,0);
	tr_gen->SetBranchStatus("pdg",1,0);
	tr_gen->SetBranchStatus("pt",1,0);
	tr_gen->SetBranchStatus("eta",1,0);

	tr_pho = (TTree*)fin->Get("multiPhotonAnalyzer/photon");
	int nPhotons;
	float pt[1000];
	float hadronicOverEm[1000];
	float phi[1000];
	int isGenMatched[1000];
	int genMomId[1000];
	tr_pho->SetBranchAddress("nPhotons",&nPhotons);
	tr_pho->SetBranchAddress("pt",pt);
	tr_pho->SetBranchAddress("hadronicOverEm",hadronicOverEm);
	tr_pho->SetBranchAddress("phi",phi);
	tr_pho->SetBranchAddress("isGenMatched",isGenMatched);
	tr_pho->SetBranchAddress("genMomId",genMomId);
	tr_pho->SetBranchStatus("*",0,0);
	tr_pho->SetBranchStatus("nPhotons",1,0);
	tr_pho->SetBranchStatus("pt",1,0);
	tr_pho->SetBranchStatus("hadronicOverEm",1,0);
	tr_pho->SetBranchStatus("phi",1,0);
	tr_pho->SetBranchStatus("isGenMatched",1,0);
	tr_pho->SetBranchStatus("genMomId",1,0);

	for(Int_t idir=0;idir<ndir;idir++)
	{
		cout <<"idir =" << idir <<" JetName ="<< DirName[idir].c_str() <<endl ;
		tr_in = (TTree*)fin->Get(Form("%s/t",DirName[idir].c_str()));

		float fentries = (float)tr_in->GetEntries();
		float weight = xSection/(fentries);

		cout<<" weight "<<weight<<" \t " << tr_in->GetName() << " entries : " << fentries << endl;

		//Declaration of leaves types
		int   nref;
		float pthat;
		float corrpt[1000];
		float jtpt[1000];
		float rawpt[1000];
		float jteta[1000];
		float jtphi[1000];
		float jty[1000];
		float jtpu[1000];
		float neutralSum[1000];
		float chargedSum[1000];
		float photonSum[1000];
		float refpt[1000];
		float refeta[1000];
		float refphi[1000];
		float refdphijt[1000];
		float refdrjt[1000];
		float refparton_pt[1000];
		int refparton_flavor[1000];
		int subid[1000];

		tr_in->SetBranchAddress("nref",&nref);
		tr_in->SetBranchAddress("pthat",&pthat);
		tr_in->SetBranchAddress("rawpt",rawpt);
		tr_in->SetBranchAddress("jtpt",jtpt);
		tr_in->SetBranchAddress("jteta",jteta);
		tr_in->SetBranchAddress("jty",jty);
		tr_in->SetBranchAddress("jtphi",jtphi);
		tr_in->SetBranchAddress("jtpu",jtpu);
		tr_in->SetBranchAddress("neutralSum",neutralSum);
		tr_in->SetBranchAddress("chargedSum",chargedSum);
		tr_in->SetBranchAddress("photonSum",photonSum);
		tr_in->SetBranchAddress("refpt",refpt);
		tr_in->SetBranchAddress("refphi",refphi);
		tr_in->SetBranchAddress("refeta",refeta);
		tr_in->SetBranchAddress("refdphijt",refdphijt);
		tr_in->SetBranchAddress("refdrjt",refdrjt);
		tr_in->SetBranchAddress("refparton_pt",refparton_pt);
		tr_in->SetBranchAddress("refparton_flavor",refparton_flavor);
		tr_in->SetBranchAddress("subid",subid);
		cout<<"get jet trees!!! "<<endl;

		//! Add Friends to the TTree
		tr_in->AddFriend(tr_ev);
		tr_in->AddFriend(tr_pho);
		tr_in->AddFriend(tr_skim);

		fout->mkdir(DirName[idir].c_str());
		fout->cd(DirName[idir].c_str());

		tr_out = new TTree("t","Jet  Response Analyzer");
		tr_out->SetMaxTreeSize(4200000000);

		int leadingPhoIndex[1000];
		int nPhogen;
		float leadingPhoPt;
		// Set output branch addresses.
	//	tr_out->Branch("pprimaryVertexFilter",&pprimaryVertexFilter,"pprimaryVertexFilter/I");
		tr_out->Branch("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA,"pPAcollisionEventSelectionPA/I");
		tr_out->Branch("pPAprimaryVertexFilter",&pPAprimaryVertexFilter,"pPAprimaryVertexFilter/I");
		tr_out->Branch("mult",&mult,"mult/I");
		tr_out->Branch("pdg",pdg,"pdg[mult]/I");
		tr_out->Branch("ptgen",ptgen,"ptgen[mult]/F");
		tr_out->Branch("etagen",etagen,"etagen[mult]/F");
		tr_out->Branch("nPhogen",&nPhogen,"nPhogen/I");
		tr_out->Branch("leadingPhoIndex",leadingPhoIndex,"leadingPhoIndex[mult]/I");
		tr_out->Branch("leadingPhoPt",&leadingPhoPt,"leadingPhoPt/F");
		tr_out->Branch("hiHF",&hiHF,"hiHF/F");
		tr_out->Branch("hiBin",&hiBin,"hiBin/I");
		tr_out->Branch("nPhotons",&nPhotons,"nPhotons/I");
		tr_out->Branch("pt",pt,"pt[nPhotons]/F");
		tr_out->Branch("hiBin",&hiBin,"hiBin/I");
		tr_out->Branch("nref",&nref,"nref/I");
		tr_out->Branch("pthat",&pthat,"pthat/F");
		tr_out->Branch("weight",&weight,"weight/F");
		tr_out->Branch("jtpt",jtpt,"jtpt[nref]/F");
		tr_out->Branch("rawpt",rawpt,"rawpt[nref]/F");
		//tr_out->Branch("corrpt",corrpt,"corrpt[nref]/F");
		tr_out->Branch("jtpu",jtpu,"jtpu[nref]/F");
		tr_out->Branch("jteta",jteta,"jteta[nref]/F");
		tr_out->Branch("jty",jty,"jty[nref]/F");
		tr_out->Branch("jtphi",jtphi,"jtphi[nref]/F");
		tr_out->Branch("neutralSum",neutralSum,"neutralSum[nref]/F");
		tr_out->Branch("chargedSum",chargedSum,"chargedSum[nref]/F");
		tr_out->Branch("photonSum",photonSum,"photonSum[nref]/F");
		tr_out->Branch("refpt",refpt,"refpt[nref]/F");
		tr_out->Branch("refeta",refeta,"refeta[nref]/F");
		tr_out->Branch("refphi",refphi,"refphi[nref]/F");
		tr_out->Branch("refdphijt",refdphijt,"refdphijt[nref]/F");
		tr_out->Branch("refdrjt",refdrjt,"refdrjt[nref]/F");
		tr_out->Branch("refparton_pt",refparton_pt,"refparton_pt[nref]/F");
		tr_out->Branch("refparton_flavor",refparton_flavor,"refparton_flavor[nref]/I");
		tr_out->Branch("subid",subid,"subid[nref]/I");


		Long64_t nbytes = 0;
		for (Long64_t i=0; i<fentries;i++) 
		{
			nbytes += tr_in->GetEntry(i);
			tr_pho->GetEntry(i);
			tr_skim->GetEntry(i);
			tr_gen->GetEntry(i);
	
			//! gen loop 
			Float_t phoPtThr = 35; // photon threshold is 35 
			Int_t nAboveThr = 0;
			Float_t leadingPt = 0;
			Int_t leadingIndex = 0;
			
			for(Int_t igen=0; igen<mult; igen++)
			{
				if(pdg[igen]!=22) continue;
				if(TMath::Abs(etagen[igen])>3) continue;
				if(ptgen[igen]<phoPtThr) continue;
				if(ptgen[igen]>leadingPt) 
				{
					leadingPt = ptgen[igen];
					leadingIndex = igen;
				}
				nAboveThr+=1;
			}
			for(Int_t igen=0; igen<mult; igen++)
                        {
				if(igen!=leadingIndex) leadingPhoIndex[igen] = 0;
				else leadingPhoIndex[igen] = 1;
			}

			leadingPhoPt = leadingPt;
			nPhogen = nAboveThr;
			if(nAboveThr==0) continue;
		

		#if 0 
			if(pthat > maxpthat) continue; 

			//if(pprimaryVertexFilter == 0) continue;// skim selection cut
			if(pPAcollisionEventSelectionPA == 0) continue;// skim selection cut
			if(pPAprimaryVertexFilter== 0) continue;// skim selection cut

			if(nPhotons == 0) continue;// no photon event

			//cout << nAboveThr << endl;
			//! loop over photons in the event
			Float_t leadingPt = 40; //minPt is 40GeV
			Int_t leadingIndex = -1;
			for(Int_t ipho = 0; ipho<nPhotons; ++ipho)
			{
				if(hadronicOverEm[ipho]>=0.1) continue;
				if(TMath::Abs(genMomId[ipho])>22) continue;
				if(isGenMatched[ipho]!=1) continue;
				if(refpt[ipho]<0) continue;
				if(pt[ipho] > leadingPt)
				{
					leadingPt = pt[ipho];
					leadingIndex = ipho;
				}
			}

			if(leadingIndex == -1) continue; // no photons above minPt

			//! jet loop
			for(int iref=0;iref<nref;iref++){
				Double_t dphi = TMath::ACos(TMath::Cos(phi[leadingIndex]-jtphi[iref]));
				if(dphi < TMath::Pi()/2){
					//        corrpt[iref]=2000.0;
					jtpt[iref]=2000.0;
					rawpt[iref]=2000.0;
					refpt[iref]=2000.0;
				} 
			}// jet delta phi loop

			#endif
			tr_out->Fill();
		}//events loop
		tr_out->Write();
	}//algorithm loop

	cout <<"finish the code!!"<<endl ;
	fout->Close();
	cout<<"XXXX"<<endl;
	//  return 0;
}

float GetXsec(float maxpt){
	float effxsec=0;
	for(int i=0; i<11; i++){
		if(TMath::Abs(maxpt - xsec[i][2]) < 1e-08){
			//effxsec = xsec[i][0] - xsec[i+1][0];
			effxsec = xsec[i][0];
			return effxsec;
		}
	}
	return  1;
}

