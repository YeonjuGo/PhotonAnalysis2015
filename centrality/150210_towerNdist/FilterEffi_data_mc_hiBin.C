// Author Yeonju Go
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "stdio.h"
#include "../HiForestAnalysis/hiForest.h"

void FilterEffi_data_mc_hiBin()
{
	
	const TCut trgCut = "HLT_HIMinBiasHfOrBSC_v1==1";
	const TCut collCut = "pcollisionEventSelection==1";
	const TCut vtxCut = "pprimaryVertexFilter==1";
	const TCut pixShapeCut = "phltPixelClusterShapeFilter";
	const TCut hf1Cut = "phfCoincFilter3==1";
	const TCut hf3Cut = "phfCoincFilter==1";
	
	const TCut runCut = "run==181611";
	const TCut lumiCut = "lumi>=1 && lumi<=895";
	const TCut eventCut = runCut && lumiCut;

	const int Ncut = 3;
	TCut totCut_data[Ncut];
        TCut totCut_mc[Ncut];
        totCut_data[0] = trgCut;
        totCut_data[1] = hf1Cut && trgCut;
        totCut_data[2] = hf3Cut && trgCut;
        totCut_mc[0] = "";
        totCut_mc[1] = hf1Cut;
        totCut_mc[2] = hf3Cut;

//	TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);

	//TFile *dataf = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/centralityDATA/merging-forest/HiForest_100_1_pKq.root");//in KNU server
	TFile *dataf = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/centralityDATA/merging-forest/HiForest_HIMinBiasUPC_Run2011_53X_run181611_09Feb2015_byYJ.root");
	//TFile *dataf = new TFile("/home/goyeonju/CMS/Files/centrality/HiForest_PbPb_minbias_DATA_20141011_53X_byKisoo.root");
	TTree *datat_evt = (TTree*) dataf -> Get("hiEvtAnalyzer/HiTree");
	TTree *datat_skim = (TTree*) dataf -> Get("skimanalysis/HltTree");
	TTree *datat_hlt = (TTree*) dataf -> Get("hltanalysis/HltTree");
	datat_evt -> AddFriend(datat_hlt);
	datat_evt -> AddFriend(datat_skim);
	double Nevt_datat = datat_evt -> GetEntries();
	cout << "# of DATA events = " << Nevt_datat << endl;

	//TFile *mcf = new TFile("/u/user/goyeonju/files/centrality/HiForest_HydjetMB_730_53XBS_merged.root");//in KNU server
	//    TFile *mcf = new TFile("/home/goyeonju/CMS/Files/centrality/HiForest_HydjetMB_730_53XBS_merged.root"); // in Korea University server
	TFile *mcf = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/centralityMC/merging_forest/centrality_PbPb_minbias_MC.root");
	TTree *mct_evt = (TTree*) mcf -> Get("hiEvtAnalyzer/HiTree");
	TTree *mct_skim = (TTree*) mcf -> Get("skimanalysis/HltTree");
	TTree *mct_hlt = (TTree*) mcf -> Get("hltanalysis/HltTree");
	mct_evt -> AddFriend(mct_hlt);
	mct_evt -> AddFriend(mct_skim);
	int Nevt_mct = mct_evt -> GetEntries();
	cout << "# of MC events = " << Nevt_mct << endl;


	//======================================
	// hiBin!!
	//======================================

	const double hiBin_bins[] = {0,2,5,10,15,20,30,40,50,60,70,80,90,100,200,300,400,500,1000,2000,3000,4000,5000};
	const int n_hiBin_bins = sizeof(hiBin_bins)/sizeof(double) - 1;

	TLine* t1 = new TLine(0,1,1000,1);
	t1->SetLineWidth(1);
	t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
	t1->SetLineColor(1); // 2 is red

	TH1D *hiBin_data[Ncut];
	TH1D *hiBin_mc[Ncut];
	for(int i=0; i<Ncut; i++)
	{
		// hiBin_data[i] = new TH1D(Form("hiBin_data%d",i), ";hiBin;Normalized Events",n_hiBin_bins, hiBin_bins); // when you look at efficiency
		hiBin_data[i] = new TH1D(Form("hiBin_data%d",i), ";hiBin;Normalized Events",210,0,210); // when you look at the whole distribution
		hiBin_data[i] -> SetMarkerStyle(20+i);
		hiBin_data[i] -> SetMarkerSize(0.7);
		hiBin_data[i] -> SetMarkerColor(kRed+i);
		hiBin_data[i] -> SetLabelSize(0.03);

		hiBin_mc[i] = (TH1D*) hiBin_data[i] -> Clone(Form("hiBin_mc%d",i));
		hiBin_mc[i] -> SetMarkerStyle(1); // marker point
		hiBin_mc[i] -> SetLineColor(2+i); //line color
		hiBin_mc[i] -> SetLabelSize(0.03);
	}
	TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
	c_temp -> cd();

	for(int i=0; i<Ncut; i++){
	
	datat_evt -> Draw("hiBin >>+ hiBin_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
	hiBin_data[0] = (TH1D*)gDirectory->Get("hiBin_data0");
	datat_evt -> Draw("hiBin >>+ hiBin_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
	hiBin_data[1] = (TH1D*)gDirectory->Get("hiBin_data1");
	datat_evt -> Draw("hiBin >>+ hiBin_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
	hiBin_data[2] = (TH1D*)gDirectory->Get("hiBin_data2");
	datat_evt -> Draw("hiBin >>+ hiBin_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
	hiBin_data[3] = (TH1D*)gDirectory->Get("hiBin_data3");
	datat_evt -> Draw("hiBin >>+ hiBin_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
	hiBin_data[4] = (TH1D*)gDirectory->Get("hiBin_data4");

	mct_evt -> Draw("hiBin >>+ hiBin_mc0");
	hiBin_mc[0] = (TH1D*)gDirectory->Get("hiBin_mc0");
	mct_evt -> Draw("hiBin >>+ hiBin_mc1","pprimaryVertexFilter==1");
	hiBin_mc[1] = (TH1D*)gDirectory->Get("hiBin_mc1");
	mct_evt -> Draw("hiBin >>+ hiBin_mc2","phltPixelClusterShapeFilter==1");
	hiBin_mc[2] = (TH1D*)gDirectory->Get("hiBin_mc2");
	mct_evt -> Draw("hiBin >>+ hiBin_mc3","phfCoincFilter3==1");
	hiBin_mc[3] = (TH1D*)gDirectory->Get("hiBin_mc3");
	mct_evt -> Draw("hiBin >>+ hiBin_mc4","pcollisionEventSelection==1");
	hiBin_mc[4] = (TH1D*)gDirectory->Get("hiBin_mc4");

	cout << "hiBin_mc[0] () Entries : " << hiBin_mc[0]->GetEntries() << endl;
	cout << "hiBin_mc[1] (pprimaryVertexFilter==1) Entries : " << hiBin_mc[1]->GetEntries() << endl;
	cout << "hiBin_mc[2] (phltPixelClusterShapeFilter==1) Entries : " << hiBin_mc[2]->GetEntries() << endl;
	cout << "hiBin_mc[3] (phfCoincFilter3==1) Entries : " << hiBin_mc[3]->GetEntries() << endl;
	cout << "hiBin_mc[4] (pcollisionEventSelection==1) Entries : " << hiBin_mc[4]->GetEntries() << endl;

	cout << "hiBin_data[0] () Entries : " << hiBin_data[0]->GetEntries() << endl;
	cout << "hiBin_data[1] (pprimaryVertexFilter==1) Entries : " << hiBin_data[1]->GetEntries() << endl;
	cout << "hiBin_data[2] (phltPixelClusterShapeFilter==1) Entries : " << hiBin_data[2]->GetEntries() << endl;
	cout << "hiBin_data[3] (phfCoincFilter3==1) Entries : " << hiBin_data[3]->GetEntries() << endl;
	cout << "hiBin_data[4] (pcollisionEventSelection==1) Entries : " << hiBin_data[4]->GetEntries() << endl;

	TLegend* l1 = new TLegend(0.3, 0.65, 0.6, 0.80, "PbPb Minbias DATA");
	// l1 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
	l1 -> AddEntry(hiBin_data[0], "No filters");
	l1 -> AddEntry(hiBin_data[1], "primay vertex filter");
	l1 -> AddEntry(hiBin_data[2], "pixel cluster shape filter");
	l1 -> AddEntry(hiBin_data[3], "HF coinc. 3 filter");
	l1 -> AddEntry(hiBin_data[4], "collision event filter");

	TLegend* l1_mc = new TLegend(0.6, 0.65, 0.9, 0.80, "PbPb Minbias MC");
	//l1 -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
	l1_mc -> AddEntry(hiBin_mc[0], "No filters");
	l1_mc -> AddEntry(hiBin_mc[1], "primay vertex filter");
	l1_mc -> AddEntry(hiBin_mc[2], "pixel cluster shape filter");
	l1_mc -> AddEntry(hiBin_mc[3], "HF coinc. 3 filter");
	l1_mc -> AddEntry(hiBin_mc[4], "collision event filter");

	TLegend* l2 = new TLegend(0.4, 0.65, 0.85, 0.80,"PbPb Minbias DATA");
	//l2 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
	l2 -> AddEntry(hiBin_data[1], "primay vertex filter");
	l2 -> AddEntry(hiBin_data[2], "pixel cluster shape filter");
	l2 -> AddEntry(hiBin_data[3], "HF coinc. 3 filter");
	l2 -> AddEntry(hiBin_data[4], "collision event filter");

	TLegend* l2_mc = new TLegend(0.4, 0.65, 0.85, 0.80, "PbPb Minbias MC");
	//l2_mc -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
	l2_mc -> AddEntry(hiBin_mc[1], "primay vertex filter");
	l2_mc -> AddEntry(hiBin_mc[2], "pixel cluster shape filter");
	l2_mc -> AddEntry(hiBin_mc[3], "HF coinc. 3 filter");
	l2_mc -> AddEntry(hiBin_mc[4], "collision event filter");

	TCanvas *c_hiBin = new TCanvas("c_hiBin", "c_hiBin", 400,400);
	c_hiBin -> SetLogy();

	double norm_data = hiBin_data[0]->Integral("width");
	double norm_mc = hiBin_mc[0]->Integral("width");
	for(int i=0; i<Ncut; i++){
		if(i==0) hiBin_data[i] -> Draw("ep");
		else hiBin_data[i] -> Draw("same&&ep");
		hiBin_mc[i] -> Draw("same&&ehist");

		//hiBin_data[i] -> Scale(1./Nevt_datat);
		//hiBin_mc[i] -> Scale(1./Nevt_mct);
		hiBin_data[i] -> Scale(1./norm_data);
		hiBin_mc[i] -> Scale(1./norm_mc);

		cout << "The integral of hiBin data " << i << " : " << hiBin_data[i]->Integral("width") << endl;
		cout << "The integral of hiBin mc " << i << " : " << hiBin_mc[i]->Integral("width") << endl;
	}

	l1 -> Draw();
	l1_mc -> Draw();

	c_hiBin -> SaveAs("pdf/hiBin.pdf");


	TH1D *hiBin_data_effi[5];
	for(int i=0; i<Ncut; i++){
		hiBin_data_effi[i] = (TH1D*)hiBin_data[i]->Clone(Form("hiBin_data_effi%d",i));
		hiBin_data_effi[i] -> SetTitle(";hiBin;Filter Efficiency");
		if(i!=0)
			hiBin_data_effi[i] -> Divide(hiBin_data[i],hiBin_data[0]);

		hiBin_data_effi[i] -> SetMarkerStyle(20+i);
		hiBin_data_effi[i] -> SetMarkerSize(0.7);
		hiBin_data_effi[i] -> SetMarkerColor(kRed+i);
		hiBin_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		//hiBin_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
	}
	TCanvas *c_hiBin_data_effi = new TCanvas("c_hiBin_data_effi", "c_hiBin_data_effi", 400, 400);
	hiBin_data_effi[1] -> Draw("elp"); 
	hiBin_data_effi[2] -> Draw("same ep"); 
	hiBin_data_effi[3] -> Draw("same ep"); 
	hiBin_data_effi[4] -> Draw("same ep"); 
	t1 -> Draw();
	l2 -> Draw();

	c_hiBin_data_effi -> SetLogx();
	c_hiBin_data_effi -> SaveAs("pdf/hiBin_data_effi.pdf");

	TH1D *hiBin_mc_effi[5];
	for(int i=0; i<Ncut; i++){
		hiBin_mc_effi[i] = (TH1D*)hiBin_mc[i]->Clone(Form("hiBin_mc_effi%d",i));
		hiBin_mc_effi[i] -> SetTitle(";hiBin;Filter Efficiency");
		if(i!=0)
			hiBin_mc_effi[i] -> Divide(hiBin_mc[i],hiBin_mc[0]);
		hiBin_mc_effi[i] -> SetMarkerStyle(20+i);
		hiBin_mc_effi[i] -> SetMarkerSize(0.7);
		hiBin_mc_effi[i] -> SetMarkerColor(2+i);
		hiBin_mc_effi[i] -> SetAxisRange(0.0,1.1,"Y");
		//hiBin_mc_effi[i] -> SetAxisRange(0.0,300.0,"X");
	}

	TCanvas *c_hiBin_mc_effi = new TCanvas("c_hiBin_mc_effi", "c_hiBin_mc_effi", 400, 400);
	hiBin_mc_effi[1] -> Draw("elp"); 
	hiBin_mc_effi[2] -> Draw("same ep"); 
	hiBin_mc_effi[3] -> Draw("same ep"); 
	hiBin_mc_effi[4] -> Draw("same ep"); 
	t1 -> Draw();
	l2_mc -> Draw();

	c_hiBin_mc_effi -> SetLogx();
	c_hiBin_mc_effi -> SaveAs("pdf/hiBin_mc_effi.pdf");

	TCanvas *c_hiBin_effi = new TCanvas("c_hiBin_effi","c_hiBin_effi", 400,400);
	c_hiBin_effi -> Divide(1,2,0);
	c_hiBin_effi -> SetLogx(); 
	c_hiBin_effi -> cd(1);
	hiBin_data_effi[0] -> Draw("ep");
	hiBin_data_effi[1] -> Draw("same ep");
	hiBin_data_effi[2] -> Draw("same ep");
	hiBin_data_effi[3] -> Draw("same ep");
	hiBin_data_effi[4] -> Draw("same ep");

	c_hiBin_effi -> cd(2);
	hiBin_mc_effi[0] -> Draw("ep");
	hiBin_mc_effi[1] -> Draw("same ep");
	hiBin_mc_effi[2] -> Draw("same ep");
	hiBin_mc_effi[3] -> Draw("same ep");
	hiBin_mc_effi[4] -> Draw("same ep");

} 
