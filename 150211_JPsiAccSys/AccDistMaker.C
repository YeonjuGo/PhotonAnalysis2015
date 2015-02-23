// modified by Yeonju
// for Acceptance systematic study

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>

const double shiftvar = -0.47; // conversion constant y=0(collision)==y=-0.47(LAB frame)  

//bool dimuCut(int, int lv_dmom0_Id, int lv_dgmom0_Id, int lv_dkid0_ch, int lv_dkid1_ch);
//bool kineCut(bool, double muPt, double muEta, double muP);
bool dimuCut(int lv_dmom0_Id, int lv_dgmom0_Id, int lv_dkid0_ch, int lv_dkid1_ch);
bool kineCut(double muPt, double muEta, double muP);
bool massCut(double lv_dimu_mass);

void AccDistMaker(char *strBinnig = "8rap9pt", bool isPrompt=true, bool isPbp=true){
    TStopwatch timer;
    timer.Start();   

	gROOT->Macro("rootlogon.C+");
	gStyle->SetCanvasDefW(800);

	char* sampleName;
	double minylab =-2.4;
	double maxylab=2.4;
	double minpt=0.0;
	double maxpt=30.0;

	//read-in tree
	TChain * ana1 = new TChain("DiAna");
	if (isPrompt){
		sampleName="PRMC_boosted";
//		ana1->Add("/home/songkyo/kyo/pPbDataSample/AcceptanceSampleKYOvtx/tot_GENonly_JPsiWithFSR_5TeV02_1st_KYOvtx_ntuple_20150126.root");
//		ana1->Add("/home/songkyo/kyo/pPbDataSample/AcceptanceSampleKYOvtx/tot_GENonly_JPsiWithFSR_5TeV02_1st_KYOvtx_ntuple_20150126_p2.root");
		ana1->Add("/home/songkyo/kyo/pPbDataSample/AcceptanceSample/pythia6_PromptJpsi_boosted_wofilters_totevtin5M_new_MuonAna_20140210.root"); //from Hyunchul
		ana1->Add("/home/songkyo/kyo/pPbDataSample/AcceptanceSample/pythia6_PromptJpsi_boosted_wofilters_totevtin5M_new_MuonAna_20140211.root"); //from Hyunchul
	}
	else {
		sampleName="NPMC_boosted";
//		ana1->Add("/home/songkyo/kyo/pPbDataSample/AcceptanceSampleKYOvtx/tot_GENonly_inclBtoJPsiMuMu_5TeV02_1st_KYOvtx_ntuple_20150127.root");
//		ana1->Add("/home/songkyo/kyo/pPbDataSample/AcceptanceSampleKYOvtx/tot_GENonly_inclBtoJPsiMuMu_5TeV02_1st_KYOvtx_ntuple_20150127_p2.root");
		ana1->Add("/home/songkyo/kyo/pPbDataSample/AcceptanceSample/pythia6_nonPromptJpsi_boosted_wofilters_totevt100M_MuonAna_20140407.root"); //from Hyunchul
	}

	const char* strName = Form("%s_%s",strBinnig,sampleName);
	std::cout << "strName: " << strName << std::endl;

	// muons (pdg +13 or -13 already)
	int dkid0_Id, dkid0_ch, dkid0_st;	
	int dkid1_Id, dkid1_ch, dkid1_st;
	double dkid0_pt, dkid0_eta, dkid0_y, dkid0_px, dkid0_py, dkid0_pz, dkid0_phi, dkid0_mass, dkid0_p;
	double dkid1_pt, dkid1_eta, dkid1_y, dkid1_px, dkid1_py, dkid1_pz, dkid1_phi, dkid1_mass, dkid1_p;
	double dimu_pt, dimu_eta, dimu_y, dimu_px, dimu_py, dimu_pz, dimu_phi, dimu_mass, dimu_p;
	// J/psi
	int dmom0_Id, dmom0_ch, dmom0_st;	
	// dgmom (not used)
	int dgmom0_Id, dgmom0_ch, dgmom0_st;
	double dmom0_pt, dmom0_eta, dmom0_y, dmom0_px, dmom0_py, dmom0_pz, dmom0_phi, dmom0_mass, dmom0_p;
	double dgmom0_pt, dgmom0_eta, dgmom0_y, dgmom0_px, dgmom0_py, dgmom0_pz, dgmom0_phi, dgmom0_mass, dgmom0_p;

	using namespace std;

	// Definition of bin
	// --- pt Bin
	Double_t ptBinsArr[] = {0.0, 3.0, 4.0, 5.0, 6.5, 7.5, 8.5, 10.0, 14.0, 30.0}; // 8rap9pt
	//Double_t ptBinsArr[] = {5.0, 6.5, 10.0, 30.0}; // 6rap3pt
	const Int_t nPtBins = sizeof(ptBinsArr)/sizeof(double)-1;
	cout << "nPtBins=" << nPtBins << endl;

	Double_t AccCentArr[] = {0.087, 0.089, 0.131, 0.196, 0.288, 0.370, 0.447, 0.539, 0.681}; // Acceptance central value

  // --- y Bin //set to 1st run (For 2nd run, will be automatically changed later)
  Double_t yBinsArr[] = {-2.4, -1.97, -1.37, -0.47, 0.43, 1.03, 1.46, 1.93, 2.4}; // 8rap9pt
  //Double_t yBinsArr[] = {-2.4, -1.97, -1.37, -0.47, 0.43, 1.03, 1.46}; // 6rap3pt
	const Int_t nYBins = sizeof(yBinsArr)/sizeof(double)-1;
	cout << "nYBins=" << nYBins << endl;

	// for 2nd run
	Double_t yBinsArr2nd[nYBins+1] = {};
	for (Int_t i=0; i<nYBins+1; i++) {
		 yBinsArr2nd[i] = -1*yBinsArr[nYBins-i];
		cout <<"yBinsArr["<<i<<"] = " <<yBinsArr[i]<<endl;
		cout <<"yBinsArr2nd["<<i<<"] = " <<yBinsArr2nd[i]<<endl;
	}
	const Int_t nYBins2nd = sizeof(yBinsArr2nd)/sizeof(double)-1;

	//read-in branches
	ana1->SetBranchAddress("dimu_eta",	&dimu_eta);
	ana1->SetBranchAddress("dimu_mass",	&dimu_mass);
	ana1->SetBranchAddress("dimu_p",	&dimu_p);
	ana1->SetBranchAddress("dimu_phi",	&dimu_phi);
	ana1->SetBranchAddress("dimu_pt",	&dimu_pt);
	ana1->SetBranchAddress("dimu_px",	&dimu_px);
	ana1->SetBranchAddress("dimu_py",	&dimu_py);
	ana1->SetBranchAddress("dimu_pz",	&dimu_pz);
	ana1->SetBranchAddress("dimu_y",	&dimu_y);

	ana1->SetBranchAddress("dgmom0_Id",	&dgmom0_Id);
	ana1->SetBranchAddress("dgmom0_ch",	&dgmom0_ch);
	ana1->SetBranchAddress("dgmom0_eta",	&dgmom0_eta);
	ana1->SetBranchAddress("dgmom0_mass",	&dgmom0_mass);
	ana1->SetBranchAddress("dgmom0_p",	&dgmom0_p);
	ana1->SetBranchAddress("dgmom0_phi",	&dgmom0_phi);
	ana1->SetBranchAddress("dgmom0_pt",	&dgmom0_pt);
	ana1->SetBranchAddress("dgmom0_px",	&dgmom0_px);
	ana1->SetBranchAddress("dgmom0_py",	&dgmom0_py);
	ana1->SetBranchAddress("dgmom0_pz",	&dgmom0_pz);
	ana1->SetBranchAddress("dgmom0_st",	&dgmom0_st);
	ana1->SetBranchAddress("dgmom0_y",	&dgmom0_y);

	ana1->SetBranchAddress("dmom0_Id",	&dmom0_Id);
	ana1->SetBranchAddress("dmom0_ch",	&dmom0_ch);
	ana1->SetBranchAddress("dmom0_eta",	&dmom0_eta);
	ana1->SetBranchAddress("dmom0_mass",	&dmom0_mass);
	ana1->SetBranchAddress("dmom0_p",	&dmom0_p);
	ana1->SetBranchAddress("dmom0_phi",	&dmom0_phi);
	ana1->SetBranchAddress("dmom0_pt",	&dmom0_pt);
	ana1->SetBranchAddress("dmom0_px",	&dmom0_px);
	ana1->SetBranchAddress("dmom0_py",	&dmom0_py);
	ana1->SetBranchAddress("dmom0_pz",	&dmom0_pz);
	ana1->SetBranchAddress("dmom0_st",	&dmom0_st);
	ana1->SetBranchAddress("dmom0_y",	&dmom0_y);

	ana1->SetBranchAddress("dkid0_Id",	&dkid0_Id);
	ana1->SetBranchAddress("dkid0_ch",	&dkid0_ch);
	ana1->SetBranchAddress("dkid0_eta",	&dkid0_eta);
	ana1->SetBranchAddress("dkid0_mass",	&dkid0_mass);
	ana1->SetBranchAddress("dkid0_p",	&dkid0_p);
	ana1->SetBranchAddress("dkid0_phi",	&dkid0_phi);
	ana1->SetBranchAddress("dkid0_pt",	&dkid0_pt);
	ana1->SetBranchAddress("dkid0_px",	&dkid0_px);
	ana1->SetBranchAddress("dkid0_py",	&dkid0_py);
	ana1->SetBranchAddress("dkid0_pz",	&dkid0_pz);
	ana1->SetBranchAddress("dkid0_st",	&dkid0_st);
	ana1->SetBranchAddress("dkid0_y",	&dkid0_y);

	ana1->SetBranchAddress("dkid1_Id",	&dkid1_Id);
	ana1->SetBranchAddress("dkid1_ch",	&dkid1_ch);
	ana1->SetBranchAddress("dkid1_eta",	&dkid1_eta);
	ana1->SetBranchAddress("dkid1_mass",	&dkid1_mass);
	ana1->SetBranchAddress("dkid1_p",	&dkid1_p);
	ana1->SetBranchAddress("dkid1_phi",	&dkid1_phi);
	ana1->SetBranchAddress("dkid1_pt",	&dkid1_pt);
	ana1->SetBranchAddress("dkid1_px",	&dkid1_px);
	ana1->SetBranchAddress("dkid1_py",	&dkid1_py);
	ana1->SetBranchAddress("dkid1_pz",	&dkid1_pz);
	ana1->SetBranchAddress("dkid1_st",	&dkid1_st);
	ana1->SetBranchAddress("dkid1_y",	&dkid1_y);

	cout << "Entries of tree : " << ana1->GetEntries() << endl;
	std::cout << "strName: " << strName << std::endl;


    ////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
    // Get Main Acceptance Value

    TString prom_st, run_st;
    if(isPrompt) prom_st = "PR"; else prom_st = "NP";
    if(isPbp) run_st = "Pbp"; else run_st = "pPb";
/*
    TFile* totAccFile = new TFile("/home/goyeonju/CMS/2015/pPbJPsiAnalysis/2015/000_fittingResult/total2Dhist_8rap9pt.root"); 
    TH2D* h2D_AccDef=(TH2D*)totAccFile->Get(Form("h2D_Acc_pt_y_%sMC_%s",prom_st.c_str(),run_st.c_str()));

    for (Int_t ipt = 0; ipt < nPtBins; ipt++) {
        if (ipt==0) {
            for (Int_t iy=1; iy< 7;iy++){
                h1D_fit_pt_y_Pbp[ipt]->SetBinContent(iy+1,0);
                h1D_fit_pt_y_Pbp[ipt]->SetBinError(iy+1,0);
                h1D_fit_pt_y_pPb[ipt]->SetBinContent(iy+1,0);
                h1D_fit_pt_y_pPb[ipt]->SetBinError(iy+1,0);
                h1D_Eff_Num_Pbp[ipt]->SetBinContent(iy+1,0);
                h1D_Eff_Num_Pbp[ipt]->SetBinError(iy+1,0);
                h1D_Eff_Num_pPb[ipt]->SetBinContent(iy+1,0);
                h1D_Eff_Num_pPb[ipt]->SetBinError(iy+1,0);
            }
        }
        else if (ipt==1 || ipt==2){
            for (Int_t iy=2; iy< 6;iy++){
                h1D_fit_pt_y_Pbp[ipt]->SetBinContent(iy+1,0);
                h1D_fit_pt_y_Pbp[ipt]->SetBinError(iy+1,0);
                h1D_fit_pt_y_pPb[ipt]->SetBinContent(iy+1,0);
                h1D_fit_pt_y_pPb[ipt]->SetBinError(iy+1,0);
                h1D_Eff_Num_Pbp[ipt]->SetBinContent(iy+1,0);
                h1D_Eff_Num_Pbp[ipt]->SetBinError(iy+1,0);
                h1D_Eff_Num_pPb[ipt]->SetBinContent(iy+1,0);
                h1D_Eff_Num_pPb[ipt]->SetBinError(iy+1,0);
            }   
        }   
        else if (ipt==3){
            for (Int_t iy=3; iy< 6;iy++){
                h1D_fit_pt_y_Pbp[ipt]->SetBinContent(iy+1,0);
                h1D_fit_pt_y_Pbp[ipt]->SetBinError(iy+1,0);
                h1D_fit_pt_y_pPb[ipt]->SetBinContent(iy+1,0);
                h1D_fit_pt_y_pPb[ipt]->SetBinError(iy+1,0);
                h1D_Eff_Num_Pbp[ipt]->SetBinContent(iy+1,0);
                h1D_Eff_Num_Pbp[ipt]->SetBinError(iy+1,0);
                h1D_Eff_Num_pPb[ipt]->SetBinContent(iy+1,0);
                h1D_Eff_Num_pPb[ipt]->SetBinError(iy+1,0);
            }   
        }   
    }   
*/

    ////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
    // Toy Acceptance Calculation
    TFile* toyFile = new TFile(Form("toyResults_pt_number100_isPrompt%d_isPbp%d.root",(int)isPrompt,(int)isPbp));
    TTree* toyTree = (TTree*)toyFile->Get("ditTree");

    double par1, par2, par3, par4; 

    toyTree->SetBranchAddress("a1",   &par1);
    toyTree->SetBranchAddress("a2",   &par2);
    toyTree->SetBranchAddress("a3",   &par3);
    toyTree->SetBranchAddress("a4",   &par4);

    TFile* outFileAcc = new TFile(Form("AccDist_isPrompt%d_isPbp%d.root",(int)isPrompt,(int)isPbp),"recreate");

    ////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	// define 1D hist
	// for Acc. Sys. study

    TH1D *h1D_Num_pt= new TH1D("h1D_Num_pt","",nPtBins,ptBinsArr);
    TH1D *h1D_Den_pt= new TH1D("h1D_Den_pt","",nPtBins,ptBinsArr);
    TH1D *h1D_Acc_pt= new TH1D("h1D_Acc_pt","",nPtBins,ptBinsArr);
    TH1D *hAccCompBin[nPtBins];
    for(int ipt=0;ipt<nPtBins;ipt++){
        //double center = hAccDef->GetBinContent(ipt+1);
        hAccCompBin[ipt] = new TH1D(Form("hAccCompBin%d",ipt),Form("Acc. Dist. of %d ptbin;Acceptance;Entries",ipt),10000, AccCentArr[ipt]-0.05, AccCentArr[ipt]+0.05);
        hAccCompBin[ipt]->Sumw2();
    }
    h1D_Num_pt->Sumw2();
    h1D_Den_pt->Sumw2();
    h1D_Acc_pt->Sumw2();

    ////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

    //// Loop Start!
    //int Ntoy = toyTree->GetEntries();
    int Ntoy = 500; 
    for(int itoy=0; itoy<Ntoy; itoy++){
        if(itoy%10==0) cout << ">>>>> itoy " << itoy << " / " << Ntoy <<  endl;
        toyTree->GetEntry(itoy);
        //cout << "par1 : " << par1 << ", par2 : " << par2 <<", par3 : " << par3 <<", par4 : " << par4 << endl;
        h1D_Num_pt-> Reset();
        h1D_Den_pt-> Reset();
        h1D_Acc_pt-> Reset();
        for(int i=0; i<ana1->GetEntries(); i++){
            // if(i%100000==0) cout << ">>>>> EVENT " << i << " / " << ana1->GetEntries() <<  endl;
            ana1->GetEntry(i);
            //weighting factor for Acc. Sys. by Yeonju
            double w_toy = 1.0;
            if(isPrompt)
                w_toy = par1/(1+TMath::Exp(par2*dimu_pt)) + par3/dimu_pt;
            else
                w_toy = par1*TMath::Erf((dimu_pt-par2)/par3) + par4;

            ////// --- cut01 : no cut
            if (!dimuCut(dmom0_Id,dgmom0_Id,dkid0_ch,dkid1_ch)) continue;

            ////// --- cut02 : for denominator
            bool yn = false;
            if (dimuCut(dmom0_Id,dgmom0_Id,dkid0_ch,dkid1_ch) && minpt<=dimu_pt && dimu_pt<maxpt && minylab<=dimu_y && dimu_y<maxylab) {yn = true;}
            if ( yn == true  ) {
                h1D_Den_pt->Fill(dimu_pt,w_toy); //yeonju
                ////// --- cut03 : for numerator
                bool yn2 = false;
                if (massCut(dimu_mass) && dimuCut(dmom0_Id,dgmom0_Id,dkid0_ch,dkid1_ch) && kineCut(dkid0_pt,dkid0_eta,dkid0_p) && kineCut(dkid1_pt,dkid1_eta,dkid1_p) && minpt<=dimu_pt && dimu_pt<maxpt && minylab<=dimu_y && dimu_y<maxylab) {yn2 = true;}
                if (yn2 == true){
                    h1D_Num_pt->Fill(dimu_pt,w_toy);//yeonju
                } // end of yn2
            } // end of yn
        } //end of loop

        // (Num/Den) to get acceptance (B : binomial error)
        h1D_Acc_pt->Divide(h1D_Num_pt,h1D_Den_pt,1,1,"B");
        for(int ipt=0;ipt<nPtBins;ipt++){
            hAccCompBin[ipt]->Fill(h1D_Acc_pt->GetBinContent(ipt+1));
           // cout << "toy number: "<< itoy << ", ptbin : " << ipt+1<< ", Acceptance : " << h1D_Acc_pt->GetBinContent(ipt+1) << endl;
        }
    } //end of itoy loop

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
	// Save the data!

    for(int ipt=0;ipt<nPtBins;ipt++){
        hAccCompBin[ipt]->Write();
    }
    outFileAcc->Close();

    timer.Stop();
    double rTime = timer.RealTime();
    cout << "real time : " << rTime << endl;
} // end of main func

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// sub-routines function 

bool dimuCut(int lv_dmom0_Id, int lv_dgmom0_Id, int lv_dkid0_ch, int lv_dkid1_ch) {
	return (TMath::Abs(lv_dmom0_Id)==443 && lv_dkid0_ch*lv_dkid1_ch<0);
}

bool kineCut(double muPt, double muEta, double muP) {
		return ( TMath::Abs(muEta) < 2.4 &&
						((TMath::Abs(muEta) < 1.3 && muPt >=3.3) ||
						 (1.3 <= TMath::Abs(muEta) && TMath::Abs(muEta) < 2.2 && muP >=2.9) ||
						 (2.2 <= TMath::Abs(muEta) && muPt >= 0.8)));
}

bool massCut(double lv_dimu_mass) {
	return ( 2.6 <= lv_dimu_mass && lv_dimu_mass < 3.5 );
}

