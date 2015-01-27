// Author : Yeonju Go
// This is modified macro from 131016_fitResolandScale/fitResolandScale.C.
// To derive residual correction fitting function.

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


void JECresidualCorr_yskim(int collision = 3, int flvOpt = 0, int genOpt = 1){
	/*  const int kHIcentral = 0; // 0-30%
	    const int kHIperipheral = 1;//30-100%
	    const int kPP = 2;
	    const int kPA = 3;
	    const int kHI010 = 4; //0-10%
	    const int kHI1030 = 5; //10-30%
	    const int kHI3050 = 6;//30-50%
	    const int kHI50100 = 7;//50-100% */

	TLegend *l1 = new TLegend(0.4365615,0.6445304,0.9577623,0.846736,NULL,"brNDC");

	TCut centCut = "";
	if ( (collision ==0) )   { // PbPb 0-30 %
		centCut = "cBin > 0 && cBin< 12";
		easyLeg(l1,"Pb+Pb 0-30%");
	}
	else if (  (collision ==1) ){  // PbPb 30-100 %
		centCut = "cBin >=12";
		easyLeg(l1,"Pb+Pb 30-100%");
	}
	else if (collision == 2 || collision == 3){ // pPb and pp 
		centCut = "";
		if (collision == 2) easyLeg(l1,"p+p");
		else easyLeg(l1, "p+Pb");
	}
	else if (  (collision ==4) ){  // PbPb 0-10 %
		centCut = "cBin > 0 && cBin < 4";
		easyLeg(l1,"Pb+Pb 0-10%");
	} else if (  (collision ==5) ){  // PbPb 10-30 %
		centCut = "cBin >= 4 && cBin < 12";
		easyLeg(l1,"Pb+Pb 10-30%");
	} else if (  (collision ==6) ){  // PbPb 30-50 %
		centCut = "cBin >= 12 && cBin < 20";
		easyLeg(l1,"Pb+Pb 30-50%");
	} else if (  (collision ==7) ){  // PbPb 50-100 %
		centCut = "cBin >= 20";
		easyLeg(l1,"Pb+Pb 50-100%");
	}    

	TH1::SetDefaultSumw2();

	// gStyle->SetOptFit(0);
	// gStyle -> SetTitleYOffset(2.35);
	gStyle -> SetOptStat(0);
	gStyle -> SetTitleYSize(0.04);


	TCut partonCut = "";
	if (flvOpt ==0 ) partonCut = "";
	else if (flvOpt == 1 ) partonCut = "refPartonFlv == 21";
	else if (flvOpt == 2 ) partonCut = "abs(refPartonFlv) < 21";
	else partonCut = "refPartonFlv < -200";

	const double ptbins[] = {30,40,50,60,70,80,90,100,120,140,160,180,200,240,280};
	const int nptbins = sizeof(ptbins)/sizeof(double) - 1;
	double AvePtBin[nptbins];

	for(int i=0;i<nptbins;i++){
		AvePtBin[i] = (ptbins[i+1]+ptbins[i])/2.0;
	}

	//###################################
	// to merge different pthat samples
	//###################################

	int nJetmax = 1000;
	float refPt[nJetmax], pt[nJetmax], eta[nJetmax], dphi[nJetmax];
	int nJet, cBin, refPartonFlv[nJetmax];
	EvtSel evtImb;
	TBranch *b_evt;

	multiTreeUtil* yJet = new multiTreeUtil();
	// PbPb
	if(collision ==1 || collision == 0  || collision == 4 || collision == 5 || collision == 6|| collision == 7){
		yJet -> addFile("/home/goyeonju/recent2013/jetAnalysis/files/yskimfiles/yskim_PbPb_pythiaHYDJET_forest_AllQCDPhotons30.root", "yJet", "");
		yJet -> addFile("/home/goyeonju/recent2013/jetAnalysis/files/yskimfiles/yskim_PbPb_pythiaHYDJET_forest_AllQCDPhotons50.root", "yJet", "");  
		yJet -> addFile("/home/goyeonju/recent2013/jetAnalysis/files/yskimfiles/yskim_PbPb_pythiaHYDJET_forest_AllQCDPhotons80.root", "yJet", "");
		yJet -> addFile("/home/goyeonju/recent2013/jetAnalysis/files/yskimfiles/yskim_PbPb_pythiaHYDJET_forest_AllQCDPhotons120.root", "yJet", "");
		yJet -> addFile("/home/goyeonju/recent2013/jetAnalysis/files/yskimfiles/yskim_PbPb_pythiaHYDJET_forest_AllQCDPhotons170.root", "yJet", "");

	} else if (collision ==3){	
#if 0
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton30to50_ak3PF.root", "yJet","",1.0 );
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton50to80_ak3PF.root", "yJet","", 1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton80to120_ak3PF.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton120to170_ak3PF.root", "yJet","",1.0 );
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton170to9999_ak3PF.root", "yJet","",1.0 );
#endif
#if 0
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton30to50_ak3PF_doJetResCorr_NoSmearing.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton50to80_ak3PF_doJetResCorr_NoSmearing.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton80to120_ak3PF_doJetResCorr_NoSmearing.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton120to170_ak3PF_doJetResCorr_NoSmearing.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton170to9999_ak3PF_doJetResCorr_NoSmearing.root", "yJet","",1.0);
#endif
#if 1
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton30to50_ak3PF_doJetResCorr_DoSmearing.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton50to80_ak3PF_doJetResCorr_DoSmearing.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton80to120_ak3PF_doJetResCorr_DoSmearing.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton120to170_ak3PF_doJetResCorr_DoSmearing.root", "yJet","",1.0);
		yJet->addFile("/u/user/goyeonju/files/yskimfiles/pA/yskim_HiForest_pPb_MIX_AllQCDPhoton170to9999_ak3PF_doJetResCorr_DoSmearing.root", "yJet","",1.0);
#endif
	} else if (collision ==2){
		//pp 
		yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_merged_allQCDPhoton30to50_genPhotonPtCut30_CMSSW538HIp2.root", "yJet", "",29329. / 29329.);
		yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_merged_allQCDPhoton50to80_genPhotonPtCut30_CMSSW538HIp2.root", "yJet", "",8098. / 87988.);
		yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_merged_allQCDPhoton80to120_genPhotonPtCut30_CMSSW538HIp2.root", "yJet", "",1680. / 96756.);
		yJet->addFile("/home/jazzitup/forestFiles/yskimmedFiles/yskim_merged_allQCDPhoton120to9999_genPhotonPtCut30_CMSSW538HIp2.root", "yJet", "", 438.   / 90972.);
	}

	yJet->AddFriend("tgj");
	// yJet->AddFriend("yPhotonTree");

	//###################################
	// flavor 
	//###################################
#if 0
	TCanvas* c_flavor = new TCanvas("c_flavor","c_flavor",400,800);
	c_flavor->Divide(1,2);
	c_flavor->cd(1);
	TH1D* hpt1 = new TH1D("hpt1",";p_{T}^{RECO};Entries",20,0,200);
	TH1D* hpt0 = (TH1D*)hpt1->Clone("hpt0");
	TH1D* hpt2 = (TH1D*)hpt1->Clone("hpt2");
	TH1D* hpt3 = (TH1D*)hpt1->Clone("hpt3");

	yJet -> Draw2(hpt0, "refPt", Form(" photonEt>40 && genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) "),"ptHatWeight*vtxCentWeight");
	yJet -> Draw2(hpt1, "refPt", Form(" photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && refPartonFlv == 21"),"ptHatWeight*vtxCentWeight");
	yJet -> Draw2(hpt2, "refPt", Form(" photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && abs(refPartonFlv)<21 "),"ptHatWeight*vtxCentWeight");
	yJet -> Draw2(hpt3, "refPt", Form(" photonEt>40 &&  genPhotonEt> 30 && abs(genMomId)<=22 && (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && refPartonFlv < -200"),"ptHatWeight*vtxCentWeight");

	handsomeTH1(hpt0,1);
	handsomeTH1(hpt1,1);
	handsomeTH1(hpt2,2);
	handsomeTH1(hpt3,4);

	hpt0->GetYaxis()->SetTitleOffset(1.8);

	hpt0->DrawCopy("hist");
	hpt1->DrawCopy("same");
	hpt2->DrawCopy("same");
	hpt3->DrawCopy("same");
	jumSun(30,0,30,7400,2);

	c_flavor->cd(2);
	hpt1->Divide(hpt0);
	hpt2->Divide(hpt0);
	hpt3->Divide(hpt0);
	hpt1->SetAxisRange(0,1,"Y");
	hpt1->SetYTitle("Ratio");
	hpt1->DrawCopy();
	hpt2->DrawCopy("same");
	hpt3->DrawCopy("same");
	jumSun(30,0,30,1,2);
#endif

	//###################################
	// to check ptHat spectrum
	//###################################
#if 0
	TCanvas* c_ptHat = new TCanvas("c_ptHat","c_ptHat",400,800);
	TH1D* hptHat = new TH1D("hptHat",";ptHat (GeV);Entries",100,0,500);

	yJet -> Draw2(hptHat, "ptHat","ptHat>0","ptHatWeight*vtxCentWeight");
	handsomeTH1(hptHat,1);
	hptHat->DrawCopy();
#endif

	//######################################################
	// pt/refpt (Jet Energy Scale) distributions for each gen jet pt bins  
	//######################################################

	TCanvas* c_JESfullDist = new TCanvas("c_JESfullDist", "pt/refPt distribution", 1200, 900); 
	makeMultiPanelCanvas(c_JESfullDist,5,4,0.0,0.0,0.2,0.15,0.02);

	TH1D* Escale[nptbins];
	double mean[nptbins], var[nptbins], resol[nptbins], resolVar[nptbins];
	for(int i=0; i < nptbins ; i++){
		c_JESfullDist -> cd(i+1);
		Escale[i] = new TH1D(Form("Escale%d_%d",collision, i) , " ;p_{T}^{RECO}/p_{T}^{GEN};Entries", 50, 0, 2);

		if ( genOpt == 1 )  {
			yJet -> Draw2(Escale[i], "pt/refPt", centCut && partonCut && Form(" (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && (refPt >= %d && refPt < %d)", (int)ptbins[i], (int)ptbins[i+1]), "ptHatWeight*vtxCentWeight"); 
			//yJet -> Draw2(Escale[i], "pt/refPt", centCut && partonCut && Form(" (abs(eta) < 3) && (dphi > 3.141592/2) && (refPt >= %d && refPt < %d)", (int)ptbins[i], (int)ptbins[i+1]), "ptHatWeight*vtxCentWeight"); 
		} else if (genOpt == 0)  {
			yJet -> Draw2(Escale[i], "pt/refPt", centCut && partonCut && Form(" (abs(eta) < 1.6) && (dphi > 7*3.141592/8.0) && (pt >= %d && pt < %d) ", (int)ptbins[i], (int)ptbins[i+1]), "ptHatWeight*vtxCentWeight");
		}
		TF1* ff = cleverGaus(Escale[i]);
		gPad->SetLogy();
		mean[i] = ff->GetParameter(1);
		var[i] = ff->GetParError(1);
		resol[i] = ff->GetParameter(2);
		resolVar[i] = ff->GetParError(2);
		cout << "mean : "<< mean[i]<< ", meanErr : " << var[i] << ", resolmean : " << resol[i] << ", resolErr : " << resolVar[i] << endl;

#if 0
		Escale[i]->Draw();	
		float dx1;    
		// ((icent==1)||(icent==4))? dx1=0.15 : dx1=0 ;
		dx1=0;
		//if ( icent == nCentBinPa )
		//  drawText(Form("E_{T}^{HF|#eta|>4} > %dGeV, ", (int)centBinPa[icent-1]), 0.12+dx1,0.929118,1,15);//yeonju 130805
		//else
		//  drawText(Form("%dGeV < E_{T}^{HF|#eta|>4} < %dGeV, ", (int)centBinPa[icent-1], (int)centBinPa[icent]), 0.12+dx1,0.929118,1,15);

		if ( i+1 == nptbins ) // This is the last bin.
			drawText(Form("p_{T}^{GEN} > %dGeV, ", (int)ptbins[i]), 0.17+dx1,0.84,1,15);//yeonju 130823
		else
			drawText(Form("%dGeV < p_{T}^{GEN} < %dGeV, ", (int)ptbins[i], (int)ptbins[i+1]), 0.17+dx1,0.84,1,12);//yeonju 130823

		//            TLegend *l1 = new TLegend(0.6365615,0.6445304,0.9577623,0.846736,NULL,"brNDC");
		//            easyLeg(l1,"p+Pb 5.02TeV");
		//    l1->AddEntry(hxjgNorm[kPADATA][icent + kPADATA*50][j],"pPb ","p");
		//                  if ( icent==1 && j==1)   l1->Draw();

#endif
	}
	c_JESfullDist -> Update();

	TCanvas* ccc = new TCanvas("ccc", "pt/refpt 30-40GeV", 400, 400); 
	ccc -> cd();
	Escale[0] -> Draw();

	//###########################################################
	// Jet Energy Scale distributions as a function of gen jet pt
	//###########################################################

	TFile* outFile = new TFile(Form("resolutionHist_collision%d.root", collision),"RECREATE");
	outFile -> cd();

	TCanvas *c_JESJER = new TCanvas("c_JESJER", ";",147,37,930,465);
	c_JESJER->Divide(2,1);
	c_JESJER->cd(1);

	TH1D* hscale = new TH1D("hscale", ";p_{T}^{GEN} (GeV);p_{T}^{RECO}/p_{T}^{GEN}", nptbins, ptbins);
	handsomeTH1(hscale,1);
	hscale -> SetAxisRange(0.9, 1.15, "Y");
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
	cout << " Jet Energy Scale fitting = C : " << par_scale[0] << ", S : " << par_scale[1] << ", N : : " << par_scale[2] << endl;

	f1->Draw("same");

	//###########################################################
	// Jet Energy Resolution distribution as a function of gen jet pt
	//###########################################################
	c_JESJER -> cd(2);
	TH1D* hresol = new TH1D("hresol", ";p_{T}^{GEN} (GeV);#sigma(p_{T}^{RECO}/p_{T}^{GEN})", nptbins, ptbins);
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
	cout << "Jet Energy Resolution fitting = C : " << par_resol[0] << ", S : " << par_resol[1] << ", N : : " << par_resol[2] << endl;

	f3->Draw("same");
	outFile -> Write();
	c_JESfullDist -> SaveAs(Form("./pdf/f1_JESfullDist_colli%d.pdf",collision));
	c_JESJER -> SaveAs(Form("./pdf/f1_JESJER_colli%d.pdf",collision));
}
