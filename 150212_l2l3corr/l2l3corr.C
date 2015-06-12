// Author : Yeonju Go
// This is modified macro from 131016_fitResolandScale/fitResolandScale.C.
// To draw l2l3 correction factor distribution.

#include "../HiForestAnalysis/hiForest.h"
#include "../gammaJetAnalysis/CutAndBinCollection2012.h"
#include <TStyle.h>
#include <TH1D.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <iostream>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TCut.h>

double myFunc(double *x, double *par){
	return par[0]+ par[1]/sqrt(x[0]) + par[2]/x[0];
}

double myFunc2(double *x, double *par){
	return sqrt(0.0567511*0.0567511+ 0.808756*0.808756/x[0] + par[0]*par[0]/(x[0]*x[0]));
}  

double myFunc3(double *x, double *par){
	return sqrt(par[0]*par[0]+ par[1]*par[1]/x[0] + par[2]*par[2]/(x[0]*x[0]));
}


void l2l3corr(int genOpt = 1, bool useFullJetTree = 0, int collision = 3, int flvOpt = 0){
	/*  const int kHIcentral = 0; // 0-30%
	    const int kHIperipheral = 1;//30-100%
	    const int kPP = 2;
	    const int kPA = 3;
	    const int kHI010 = 4; //0-10%
	    const int kHI1030 = 5; //10-30%
	    const int kHI3050 = 6;//30-50%
	    const int kHI50100 = 7;//50-100% */

	TLegend *l1 = new TLegend(0.4365615,0.6445304,0.9577623,0.846736,NULL,"brNDC");

	TH1::SetDefaultSumw2();

	// gStyle->SetOptFit(0);
	// gStyle -> SetTitleYOffset(2.35);
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.04);

	//const double ptbins[] = {30,40,50,60,80,100,140,180,280};
	const double ptbins[] = {30,40,50,60,70,80,90,100,120,140,160,180,200,240,280};
	const int nptbins = sizeof(ptbins)/sizeof(double) - 1;
	double AvePtBin[nptbins];

	for(int i=0;i<nptbins;i++){
		AvePtBin[i] = (ptbins[i+1]+ptbins[i])/2.0;
	}

	//###################################
	// to merge different pthat samples
	//###################################

	int nJetmax = 100;
	float refPt[nJetmax], pt[nJetmax], eta[nJetmax], dphi[nJetmax];
	int nJet, cBin, refPartonFlv[nJetmax];
	EvtSel evtImb;
	TBranch *b_evt;
	TString treeName = "yJet";
	if(useFullJetTree==1) treeName = "fullJet";

	multiTreeUtil* yJet = new multiTreeUtil();
	if (collision ==3){	
#if 1
		// 2015/02/12 Thur. l2l3 correction study
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton30to50_akPu3PF_AfterResCorr_final_l2l3corrTest.root",treeName,"");
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton50to80_akPu3PF_AfterResCorr_final_l2l3corrTest.root",treeName ,"");
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton80to120_akPu3PF_AfterResCorr_final_l2l3corrTest.root",treeName ,"");
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton120to170_akPu3PF_AfterResCorr_final_l2l3corrTest.root",treeName ,"");
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton170to9999_akPu3PF_AfterResCorr_final_l2l3corrTest.root",treeName ,"");
#endif

	} 

	yJet->AddFriend("tgj");
	// yJet->AddFriend("yPhotonTree");


	//###################################
	// to check ptHat spectrum
	//###################################
#if 1
	if(useFullJetTree==0){
		TCanvas* c_ptHat = new TCanvas("c_ptHat","c_ptHat",400,400);
		TH1D* hptHat = new TH1D("hptHat",";ptHat (GeV);Entries",100,0,500);
		//yJet -> Draw2(hptHat, "ptHat","ptHat>0");
		yJet -> Draw2(hptHat, "ptHat","ptHat>0","ptHatWeight*vtxCentWeight");
		handsomeTH1(hptHat,1);
		hptHat->DrawCopy();
	}
#endif

	//######################################################
	// l2l3 correction factors as a function of pt 
	//######################################################

	TH1D* h1D_l2l3_gen[nptbins];
	TH1D* h1D_l2l3_reco[nptbins];
	TH1D* h2D_pt_reco[nptbins];
	for(int ipt=0;ipt<nptbins;ipt++){
		h1D_l2l3_gen[ipt] = new TH1D(Form("h1D_l2l3_gen%d",ipt),"",1000,0,2);
		h1D_l2l3_reco[ipt] = new TH1D(Form("h1D_l2l3_reco%d",ipt),"",1000,0,2);
	}
	
	TCanvas* c_gen_dist = new TCanvas("c_gen_dist", "l2l3 distribution for each gen pt bin", 1200, 900); 
	makeMultiPanelCanvas(c_gen_dist,5,4,0.0,0.0,0.2,0.15,0.02);

	double mean[nptbins], var[nptbins], resol[nptbins], resolVar[nptbins];
	double mean_reco[nptbins], var_reco[nptbins], resol_reco[nptbins], resolVar_reco[nptbins];
	for(int ipt=0;ipt<nptbins;ipt++){
		c_gen_dist->cd(ipt+1);
		yJet -> Draw2(h1D_l2l3_gen[ipt], "l2l3corr", Form("(abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && (refPt >= %d && refPt < %d)", (int)ptbins[ipt], (int)ptbins[ipt+1]), "ptHatWeight*vtxCentWeight"); 
		yJet -> Draw2(h1D_l2l3_reco[ipt], "l2l3corr", Form("(abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && (pt >= %d && pt < %d)", (int)ptbins[ipt], (int)ptbins[ipt+1]), "ptHatWeight*vtxCentWeight"); 	

		TF1* ff_gen = cleverGaus(h1D_l2l3_gen[ipt]);
		gPad->SetLogy();
		mean[ipt] = ff_gen->GetParameter(1);
		var[ipt] = ff_gen->GetParError(1);
		resol[ipt] = ff_gen->GetParameter(2);
		resolVar[ipt] = ff_gen->GetParError(2);
		cout << "mean : "<< mean[ipt]<< ", meanErr : " << var[ipt] << ", resolmean : " << resol[ipt] << ", resolErr : " << resolVar[ipt] << endl;

		TF1* ff_reco = cleverGaus(h1D_l2l3_reco[ipt]);
		gPad->SetLogy();
		mean_reco[ipt] = ff_reco->GetParameter(1);
		var_reco[ipt] = ff_reco->GetParError(1);
		resol_reco[ipt] = ff_reco->GetParameter(2);
		resolVar_reco[ipt] = ff_reco->GetParError(2);
		cout << "mean : "<< mean_reco[ipt]<< ", meanErr : " << var_reco[ipt] << ", resolmean : " << resol_reco[ipt] << ", resolErr : " << resolVar_reco[ipt] << endl;
	}
	c_gen_dist -> Update();

	TCanvas* ccc = new TCanvas("ccc", "pt/refpt 30-40GeV", 400, 400); 
	ccc -> cd();

	//###########################################################
	// Jet Energy Scale distributions as a function of gen jet pt
	//###########################################################
#if 0
	TFile* outFile = new TFile(Form("resolutionHist_collision%d.root", collision),"RECREATE");
	outFile -> cd();

	TCanvas *c_JESJER = new TCanvas("c_JESJER", ";",37,147,465,930);
	c_JESJER->Divide(1,2);
	c_JESJER->cd(1);

	TH1D* hscale = new TH1D("hscale", ";p_{T}^{GEN} (GeV);p_{T}^{RECO}/p_{T}^{GEN}", nptbins, ptbins);
	if(genOpt==0) hscale->SetXTitle("p_{T}^{RECO} (GeV)");
	handsomeTH1(hscale,1);
	hscale -> SetAxisRange(0.8, 1.2, "Y");
	hscale -> Draw();
	jumSun(30,1,200,1);

	TLegend *l3=new TLegend(0,0,0.4490239,0.08695652,NULL,"brNDC");
	l3->SetTextFont(42);
	l3->SetTextSize(0.04);
	l3->SetFillColor(0);
	l3->SetLineColor(0);

	for(int i=0; i < nptbins ; i++){ 
		hscale -> SetBinContent(i+1, mean[i]);
		//hscale -> SetBinError(i+1, 0.00001);
		hscale -> SetBinError(i+1, var[i]);
	}
	
	TF1 *f1 = new TF1("f1", myFunc, 30, 300, 3);
	f1 -> SetParameters(0.9,0.8,0.001);
	f1 -> SetParNames("C1", "S1", "N1");

	TF1 *f2;
	f2 = new TF1("f2", myFunc2, 30, 300, 3);
	f2 -> SetParameters(0.03, 0.8, 0.01);
	f2 -> SetParNames("C2", "S2", "N2");

	TF1 *f3;
	f3 = new TF1("f3", myFunc3, 30, 300, 3);
	f3 -> SetParameters(0.03, 0.8, 0.01);
	f3 -> SetParNames("C3", "S3", "N3");


	double par_scale[3];

	hscale -> Fit("f1", "RLL");
	f1->GetParameters(par_scale);
	cout << " Jet Energy Scale fitting = C : " << par_scale[0] << ", S : " << par_scale[1] << ", N : " << par_scale[2] << endl;

	f1->Draw("same");

	//###########################################################
	// Jet Energy Resolution distribution as a function of gen jet pt
	//###########################################################
	c_JESJER -> cd(2);
	TH1D* hresol = new TH1D("hresol", ";p_{T}^{GEN} (GeV);#sigma(p_{T}^{RECO}/p_{T}^{GEN})", nptbins, ptbins);
	if(genOpt==0) hresol->SetXTitle("p_{T}^{RECO} (GeV)");
	handsomeTH1(hresol,1);
	hresol->SetAxisRange(0.0, 0.3, "Y");
	hresol->GetYaxis()->CenterTitle();
	hresol->Draw();
	jumSun(30,0.0,30,0.3,2);

	for(int i=0; i < nptbins ; i++){ 
		hresol -> SetBinContent(i+1, resol[i]);
		//hresol -> SetBinError(i+1, 0.00001);
		hresol -> SetBinError(i+1, resolVar[i]);
	}

	double par_resol[3];

	hresol -> Fit("f3", "RLL");
	f3->GetParameters(par_resol);
	cout << "Jet Energy Resolution fitting = C : " << par_resol[0] << ", S : " << par_resol[1] << ", N : " << par_resol[2] << endl;

	f3->Draw("same");
	outFile -> Write();
	c_JESfullDist -> SaveAs(Form("./pdf/f1_JESfullDist_colli%d.pdf",collision));
	c_JESJER -> SaveAs(Form("./pdf/f1_JESJER_colli%d.pdf",collision));
#endif
}
