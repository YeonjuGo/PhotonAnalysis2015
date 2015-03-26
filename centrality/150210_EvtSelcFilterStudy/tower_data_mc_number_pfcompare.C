// Author Yeonju Go
// last modification : 2015/01/27 
// 
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
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
#include "TLatex.h"
#include "stdio.h"
#include "../../HiForestAnalysis/hiForest.h"
#include "../../gammaJetAnalysis/commonUtility.h"

void tower_data_mc_number_pfcompare(float etThr=0.5)
{
	const int Ncut = 3;
	const TCut runCut = "run==181611";
	const TCut lumiCut = "lumi>=1 && lumi<=895";
	// const TCut eventCut = "";
	const TCut eventCut = runCut && lumiCut;
//	const TCut etaCut = "abs(eta)>2.87 && abs(eta)<5.2";
//	const TCut etCut = "et>1.4";
	const TCut etaCut = "";
	const TCut etCut = "";
	const TCut trgCut = "HLT_HIMinBiasHfOrBSC_v1==1";
	const TCut hf1Cut = "phfCoincFilter==1";
	const TCut hf3Cut = "phfCoincFilter3==1";
	TCut totCut_data[Ncut];
	TCut totCut_mc[Ncut];
	totCut_data[0] = eventCut && etaCut && etCut && trgCut;
	totCut_data[1] = hf1Cut && eventCut && etaCut && etCut && trgCut;
	totCut_data[2] = hf3Cut && eventCut && etaCut && etCut && trgCut;
	totCut_mc[0] = etaCut && etCut;
	totCut_mc[1] = hf1Cut && etaCut && etCut;
	totCut_mc[2] = hf3Cut && etaCut && etCut;

	double color[Ncut] = {12,8,9}; //{,41,46};
	double marker[Ncut] = {20,22,29}; //{,33,34};
	string FilterName[Ncut] = {"no evt filter","phfCoincFilter","phfCoincFilter3"};
	//string FilterName[Ncut] = {"no evt filter","pprimaryVertexFilter","phltPixelClusterShapeFilter", "phfCoincFilter3", "pcollisionEventSelection"};

	//TH1::SetDefaultSumw2();
	gStyle -> SetOptStat(0);

	TLatex* latex = new TLatex();
	latex->SetNDC();
	latex->SetTextAlign(12);
	latex->SetTextSize(0.04);
	// ===================================================================================
	// Get Trees from data & mc files.
	// ===================================================================================
	TFile *dataf = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/centralityDATA/merging-forest/HiForest_HIMinBiasUPC_Run2011_53X_run181611_09Feb2015_byYJ.root");
	TTree *datat_evt = (TTree*) dataf -> Get("hiEvtAnalyzer/HiTree");
	TTree *datat_skim = (TTree*) dataf -> Get("skimanalysis/HltTree");
	TTree *datat_hlt = (TTree*) dataf -> Get("hltanalysis/HltTree");
	TTree *datat = (TTree*) dataf -> Get("rechitanalyzer/tower");
	TTree *datat_recTower1 = (TTree*) dataf -> Get("rechitanalyzer/hf");
	TTree *datat_recTower2 = (TTree*) dataf -> Get("rechitanalyzer/hbhe");
	TTree *datat_recTower3 = (TTree*) dataf -> Get("rechitanalyzer/ee");
	TTree *datat_recTower4 = (TTree*) dataf -> Get("rechitanalyzer/eb");
	TTree *datat_recTower5 = (TTree*) dataf -> Get("rechitanalyzer/ntEvent");
	datat -> AddFriend(datat_hlt);
	datat -> AddFriend(datat_evt);
	datat -> AddFriend(datat_skim);
	datat -> AddFriend(datat_recTower1);
	datat -> AddFriend(datat_recTower2);
	datat -> AddFriend(datat_recTower3);
	datat -> AddFriend(datat_recTower4);
	datat -> AddFriend(datat_recTower5);

	TFile *mcf = new TFile("/u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/centralityMC/merging_forest/centrality_PbPb_minbias_MC.root");
	//TFile *mcf = new TFile("/home/goyeonju/CMS/Files/centrality/centrality_PbPb_minbias_MC_53X_byYJ.root");
	//TFile *mcf = new TFile("/home/goyeonju/CMS/Files/centrality/HiForest_HydjetMB_730_53XBS_merged.root");
	TTree *mct_evt = (TTree*) mcf -> Get("hiEvtAnalyzer/HiTree");
	TTree *mct_skim = (TTree*) mcf -> Get("skimanalysis/HltTree");
	TTree *mct_hlt = (TTree*) mcf -> Get("hltanalysis/HltTree");
	TTree *mct = (TTree*) mcf -> Get("rechitanalyzer/tower");
	TTree *mct_recTower1 = (TTree*) mcf -> Get("rechitanalyzer/hf");
	TTree *mct_recTower2 = (TTree*) mcf -> Get("rechitanalyzer/hbhe");
	TTree *mct_recTower3 = (TTree*) mcf -> Get("rechitanalyzer/ee");
	TTree *mct_recTower4 = (TTree*) mcf -> Get("rechitanalyzer/eb");
	TTree *mct_recTower5 = (TTree*) mcf -> Get("rechitanalyzer/ntEvent");
	mct -> AddFriend(mct_hlt);
	mct -> AddFriend(mct_evt);
	mct -> AddFriend(mct_skim);
	mct -> AddFriend(mct_recTower1);
	mct -> AddFriend(mct_recTower2);
	mct -> AddFriend(mct_recTower3);
	mct -> AddFriend(mct_recTower4);
	mct -> AddFriend(mct_recTower5);

	// ===============================================================================================
	// [et] Define et histograms (data/mc , pf tree/rechit tree) 
	// ===============================================================================================
	cout << "LET'S DEFINE HISTOGRAMS" << endl;

	TH1F* data_rec_n[Ncut];
	TH1F* mc_rec_n[Ncut];
	TH1F* ratio_rec_n[Ncut];// ratio = data/mc

	for(int i=0; i<Ncut; i++)
	{
		data_rec_n[i] = new TH1F(Form("data_rec_n%d",i), ";# of HF towers;Events",200,0.001,1000);
		data_rec_n[i] -> SetMarkerStyle(20);
		data_rec_n[i] -> SetMarkerSize(0.4);
		data_rec_n[i] -> SetMarkerColor(color[i]);
		data_rec_n[i] -> SetLineColor(color[i]);
		data_rec_n[i] -> SetLabelSize(0.03);

		mc_rec_n[i] = (TH1F*)data_rec_n[i]->Clone(Form("mc_rec_n%d",i));
		mc_rec_n[i] -> SetMarkerStyle(marker[i]);
		mc_rec_n[i] -> SetMarkerSize(0.9);
		mc_rec_n[i] -> SetMarkerColor(color[i]); //marker color
		mc_rec_n[i] -> SetLineColor(color[i]); //line color
		mc_rec_n[i] -> SetLabelSize(0.03);

		ratio_rec_n[i] = (TH1F*)data_rec_n[i]->Clone(Form("ratio_rec_n%d",i));
		ratio_rec_n[i] -> SetTitle(";# of HF towers; (evt of filter)/(evt of no filter)"); 
		ratio_rec_n[i] -> SetMarkerStyle(20);
		ratio_rec_n[i] -> SetMarkerSize(0.5);
		ratio_rec_n[i] -> SetMarkerColor(color[i]);
		ratio_rec_n[i] -> SetLineColor(color[i]);
		ratio_rec_n[i] -> SetLabelSize(0.03);
		ratio_rec_n[i] -> GetXaxis() -> SetRangeUser(0,100);
		//ratio_rec_n[i] -> SetAxisRange(0.5,1.5,"Y");

	}

	// ===============================================================================================
	// [rechit] Fill histogram by using <Draw> function in TTree.
	// ===============================================================================================
	cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;

	for(int i=0; i<Ncut; i++){
		datat -> Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.et > %f)>> %s",etThr,data_rec_n[i]->GetName()),totCut_data[i]);
		mct -> Draw(Form("Sum$(abs(tower.eta) > 2.87 && tower.et > %f)>> %s",etThr,mc_rec_n[i]->GetName()),totCut_mc[i]);
	}

	TCanvas *c_mc = new TCanvas("c_mc","c_mc", 500,350*2);
	c_mc->Divide(1,2);
	TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
	easyLeg(leg);
	for(int i=0; i<Ncut; i++){
		int cutBinFrom = data_rec_n[i]->FindBin(100); 
                int cutBinTo = data_rec_n[i]->FindBin(800);
		data_rec_n[i] -> Scale(mc_rec_n[i]->Integral(cutBinFrom,cutBinTo)/data_rec_n[i]->Integral(cutBinFrom,cutBinTo));

		c_mc->cd(1);
		if(i==0) { data_rec_n[i]->Draw("p"); mc_rec_n[i]->Draw("same hist");}
		else {
			mc_rec_n[i]->Draw("same hist");
			ratio_rec_n[i]->Divide(mc_rec_n[i],mc_rec_n[0],1.0,1.0,"B");
			ratio_rec_n[i] -> GetYaxis() -> SetRangeUser(0.5,1.5);
		}
		cout << i << "th integral : " << mc_rec_n[i]->Integral() << endl;
		cout << i << "th integral in the range (100,800) : " << mc_rec_n[i]->Integral(cutBinFrom, cutBinTo) << endl;
	}
	c_mc->cd(2);
	ratio_rec_n[2]->Draw("hist");
	ratio_rec_n[1]->Draw("same hist");
//	leg->SetFillColor(0);
//	leg->SetTextFont(43);
//	leg->SetTextSize(15);
	leg->SetHeader("rechitTowers");
	leg->AddEntry(data_rec_n[0],"DATA","p");
	leg->AddEntry(mc_rec_n[0],"MC no filter","l");
	leg->AddEntry(mc_rec_n[1],Form("MC %s",hf1Cut.GetTitle()),"l");
	leg->AddEntry(mc_rec_n[2],Form("MC %s",hf3Cut.GetTitle()),"l");
	c_mc->cd(1);
	leg->Draw();

	TFile* outFile = new TFile("outfile.root","recreate");
	outFile->cd();
	for(int i=0;i<Ncut;i++){
		mc_rec_n[i]->Write();
		data_rec_n[i]->Write();
		ratio_rec_n[i]->Write();
	}
	c_mc->SaveAs(Form("pdf/mc_hfcoin_etThr%f.pdf",etThr));
	outFile->Close();

#if 0
	// ===============================================================================================
	// [# of HF tower] Draw histograms in Canvas. 
	// ===============================================================================================
	cout << "LET'S DRAW HISTOGRAMS IN CANVAS" << endl;
	TCanvas *c_tot = new TCanvas("c_tot","c_tot", 1500,350);
	makeMultiPanelCanvas(c_tot,5,1,0.0,0.0,0.2,0.15,0.02);
	for(int i=0; i<Ncut; i++){
		c_tot->cd(i+1);
		mc_rec_n[i]->DrawCopy("hist");
		data_rec_n[i]->Draw("same");
		if(i==0)latex->DrawLatex(0.46,0.88,"rechitTowers");
		latex->DrawLatex( 0.46, 0.78, Form("%s", FilterName[i].c_str()) );
	}

	TCanvas *c_mc = new TCanvas("c_mc", "c_mc", 400,800);
	//    makeMultiPanelCanvas(c_mc,1,2,0.0,0.0,0.2,0.15,0.02);
	c_mc -> Divide(1,2);

	c_mc->cd(1);
	gPad->SetLogy();
	mc_rec_n[0]->GetXaxis()->SetRangeUser(0.,250.);
	mc_rec_n[0] -> Draw("ep");
	mc_rec_n[1] -> Draw("same ep");
	mc_rec_n[2] -> Draw("same ep");
	mc_rec_n[3] -> Draw("same ep");
	mc_rec_n[4] -> Draw("same ep");
	//  latex->DrawLatex(0.46,0.68,"rechitTowers MC");

	c_mc->cd(2);
	gPad->SetLogy();
	TH1F* hFilterRatio = (TH1F*)mc_rec_n[0]->Clone("hFilterRatio");
	hFilterRatio->GetYaxis()->SetTitle("Filter Efficiency");
	for(int i=1;i<Ncut;i++){
		hFilterRatio->Reset();
		hFilterRatio->GetXaxis()->SetRangeUser(0.,250.);
		hFilterRatio->Divide(mc_rec_n[i],mc_rec_n[0]);
		hFilterRatio->SetMarkerColor(color[i]);
		hFilterRatio->SetMarkerStyle(marker[i]);
		if(i==1) hFilterRatio->DrawCopy();
		hFilterRatio->DrawCopy("same");
	}

#endif
} 
