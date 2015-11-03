/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */
#include "gedPhotonUtility.h"

void draw_gedPhoton_effi(float ptThr=30, condition cond_ = noC )
{
	gStyle->SetOptStat(0);
	gStyle->SetTitleXSize(0.08);
	gStyle->SetTitleYSize(0.08);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetHistLineWidth(1);
	TH1::SetDefaultSumw2();

	TString dirName = getCondDirName(cond_);
	TString outSuffix = getCondSuffix(cond_);

	const int nReco= 2;//number of files
	const char* fName1 = Form("histFiles/gedPhotonHist_effi_ggHiNtuplizer_ptThr%d%s.root",(int)ptThr,outSuffix.Data());
	const char* fName2 = Form("histFiles/gedPhotonHist_effi_ggHiNtuplizerGED_ptThr%d%s.root",(int)ptThr,outSuffix.Data());

	TFile* fin[nReco];
	fin[0] = new TFile(fName1, "READ");
	fin[1] = new TFile(fName2, "READ");
	std::cout << "input HiForest 1 : " << fin[0]->GetName() << std::endl;
	std::cout << "input HiForest 2 : " << fin[1]->GetName() << std::endl;

	// histos[4][5][5];
	// [4][][] : eta cuts : no cut, Barrel, Endcap, Endcap2
	// [][5][] : hiBin cuts : all centralities, 0-10%, 10-30%, 30-50%, 50-100%

	TH1D* h_detIso_effi[nEtaCut][nRadius][2]; //as a fnc. of hiBin
	TH1D* h_detIso_frac[nEtaCut][nRadius][2]; //as a fnc. of hiBin
	TH1D* h_detIso_signal[nEtaCut][nRadius][2]; //as a fnc. of hiBin
	TH1D* h_pfIso_effi[nEtaCut][nRadius][2]; //as a fnc. of hiBin
	TH1D* h_pfIso_frac[nEtaCut][nRadius][2]; //as a fnc. of hiBin
	TH1D* h_pfIso_signal[nEtaCut][nRadius][2]; //as a fnc. of hiBin
	TH1D* h_pfVsIso_effi[nEtaCut][nRadius][2]; //as a fnc. of hiBin
	TH1D* h_pfVsIso_frac[nEtaCut][nRadius][2]; //as a fnc. of hiBin
	TH1D* h_pfVsIso_signal[nEtaCut][nRadius][2]; //as a fnc. of hiBin

	TCanvas* c_effi_frac[nEtaCut][nRadius];
	TCanvas* c_effi_frac_etaTogether[nRadius];
	int histColor[2] = {1,1};
	TLegend* leg = new TLegend(0.6,0.2,0.8,0.6);
	legStyle(leg);
	for(int iR=2; iR<nRadius; ++iR){
		c_effi_frac_etaTogether[iR] = new TCanvas(Form("c_effi_frac_etaTogether_R%d",iR),"", 1000,600);
		c_effi_frac_etaTogether[iR]->Divide(3,3);
		for(int j=1; j<nEtaCut; ++j){
			c_effi_frac[j][iR] = new TCanvas(Form("c_effi_frac_R%d_eta%d",iR,j),"", 400,600);
			c_effi_frac[j][iR]->Divide(1,3);

			h_detIso_effi[j][iR][1]->SetMarkerStyle(20);
			h_detIso_frac[j][iR][1]->SetMarkerStyle(20);
			h_detIso_signal[j][iR][1]->SetMarkerStyle(20);
			h_pfIso_effi[j][iR][1]->SetMarkerStyle(20);
			h_pfIso_frac[j][iR][1]->SetMarkerStyle(20);
			h_pfIso_signal[j][iR][1]->SetMarkerStyle(20);
			h_pfVsIso_effi[j][iR][1]->SetMarkerStyle(20);
			h_pfVsIso_frac[j][iR][1]->SetMarkerStyle(20);
			h_pfVsIso_signal[j][iR][1]->SetMarkerStyle(20);

			for(int ifile=0; ifile<nReco; ++ifile){
				TString name = Form("R%d_eta%d_momId0",iR,j);
				h_detIso_effi[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_detIso_effi_%s",name.Data()));
				h_detIso_frac[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_detIso_frac_%s",name.Data()));
				h_detIso_signal[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_detIso_signal_%s",name.Data()));
				h_pfIso_effi[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_pfIso_effi_%s",name.Data()));
				h_pfIso_frac[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_pfIso_frac_%s",name.Data()));
				h_pfIso_signal[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_pfIso_signal_%s",name.Data()));
				h_pfVsIso_effi[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_pfVsIso_effi_%s",name.Data()));
				h_pfVsIso_frac[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_pfVsIso_frac_%s",name.Data()));
				h_pfVsIso_signal[j][iR][ifile] = (TH1D*) fin[ifile]->Get(Form("h_pfVsIso_signal_%s",name.Data()));
		
				h_detIso_effi[j][iR][ifile]->SetMarkerColor(kOrange+histColor[ifile]);
				h_detIso_frac[j][iR][ifile]->SetMarkerColor(kOrange+histColor[ifile]);
				h_detIso_signal[j][iR][ifile]->SetMarkerColor(kOrange+histColor[ifile]);
				h_pfIso_effi[j][iR][ifile]->SetMarkerColor(kViolet+histColor[ifile]);
				h_pfIso_frac[j][iR][ifile]->SetMarkerColor(kViolet+histColor[ifile]);
				h_pfIso_signal[j][iR][ifile]->SetMarkerColor(kViolet+histColor[ifile]);
				h_pfVsIso_effi[j][iR][ifile]->SetMarkerColor(kGreen+histColor[ifile]);
				h_pfVsIso_frac[j][iR][ifile]->SetMarkerColor(kGreen+histColor[ifile]);
				h_pfVsIso_signal[j][iR][ifile]->SetMarkerColor(kGreen+histColor[ifile]);
				
				h_detIso_effi[j][iR][ifile]->SetLineColor(kOrange+histColor[ifile]);
				h_detIso_frac[j][iR][ifile]->SetLineColor(kOrange+histColor[ifile]);
				h_detIso_signal[j][iR][ifile]->SetLineColor(kOrange+histColor[ifile]);
				h_pfIso_effi[j][iR][ifile]->SetLineColor(kViolet+histColor[ifile]);
				h_pfIso_frac[j][iR][ifile]->SetLineColor(kViolet+histColor[ifile]);
				h_pfIso_signal[j][iR][ifile]->SetLineColor(kViolet+histColor[ifile]);
				h_pfVsIso_effi[j][iR][ifile]->SetLineColor(kGreen+histColor[ifile]);
				h_pfVsIso_frac[j][iR][ifile]->SetLineColor(kGreen+histColor[ifile]);
				h_pfVsIso_signal[j][iR][ifile]->SetLineColor(kGreen+histColor[ifile]);
			
				h_detIso_effi[j][iR][ifile]->GetYaxis()->SetRangeUser(0.5,1);
				h_detIso_frac[j][iR][ifile]->GetYaxis()->SetRangeUser(0.5,1);
				h_detIso_signal[j][iR][ifile]->GetYaxis()->SetRangeUser(0.3,1);
			}//ifile
			if(iR==2 && j==1) {
				leg->AddEntry(h_detIso_effi[j][iR][0], "OLD detIso");
				leg->AddEntry(h_detIso_effi[j][iR][1], "GED detIso");
				leg->AddEntry(h_pfIso_effi[j][iR][0], "OLD pfIso");
				leg->AddEntry(h_pfIso_effi[j][iR][1], "GED pfIso");
				leg->AddEntry(h_pfVsIso_effi[j][iR][0], "OLD pfVsIso");
				leg->AddEntry(h_pfVsIso_effi[j][iR][1], "GED pfVsIso");
			}
			// seperate canvas for each eta
			c_effi_frac[j][iR]->cd(1);
			h_detIso_effi[j][iR][0]->Draw("hist");
			h_detIso_effi[j][iR][1]->Draw("pl same");
			h_pfIso_effi[j][iR][0]->Draw("hist same");
			h_pfIso_effi[j][iR][1]->Draw("pl same");
			h_pfVsIso_effi[j][iR][0]->Draw("hist same");
			h_pfVsIso_effi[j][iR][1]->Draw("hist same");
			leg->Draw("same");
			c_effi_frac[j][iR]->cd(2);
			h_detIso_frac[j][iR][0]->Draw("hist");
			h_detIso_frac[j][iR][1]->Draw("hist same");
			h_pfIso_frac[j][iR][0]->Draw("hist same");
			h_pfIso_frac[j][iR][1]->Draw("hist same");
			h_pfVsIso_frac[j][iR][0]->Draw("hist same");
			h_pfVsIso_frac[j][iR][1]->Draw("hist same");

			c_effi_frac[j][iR]->cd(3);
			h_detIso_signal[j][iR][0]->Draw("hist");
			h_detIso_signal[j][iR][1]->Draw("hist same");
			h_pfIso_signal[j][iR][0]->Draw("hist same");
			h_pfIso_signal[j][iR][1]->Draw("hist same");
			h_pfVsIso_signal[j][iR][0]->Draw("hist same");
			h_pfVsIso_signal[j][iR][1]->Draw("hist same");


			// combined canvas in terms of eta
			c_effi_frac_etaTogether[iR]->cd(j);
			h_detIso_effi[j][iR][0]->Draw("hist");
			h_detIso_effi[j][iR][1]->Draw("pl same");
			h_pfIso_effi[j][iR][0]->Draw("hist same");
			h_pfIso_effi[j][iR][1]->Draw("pl same");
			h_pfVsIso_effi[j][iR][0]->Draw("hist same");
			h_pfVsIso_effi[j][iR][1]->Draw("pl same");
			if(j==1) leg->Draw("same");
			c_effi_frac_etaTogether[iR]->cd(j+3);
			h_detIso_frac[j][iR][0]->Draw("hist");
			h_detIso_frac[j][iR][1]->Draw("pl same");
			h_pfIso_frac[j][iR][0]->Draw("hist same");
			h_pfIso_frac[j][iR][1]->Draw("pl same");
			h_pfVsIso_frac[j][iR][0]->Draw("hist same");
			h_pfVsIso_frac[j][iR][1]->Draw("pl same");

			c_effi_frac_etaTogether[iR]->cd(j+6);
			h_detIso_signal[j][iR][0]->Draw("hist");
			h_detIso_signal[j][iR][1]->Draw("pl same");
			h_pfIso_signal[j][iR][0]->Draw("hist same");
			h_pfIso_signal[j][iR][1]->Draw("pl same");
			h_pfVsIso_signal[j][iR][0]->Draw("hist same");
			h_pfVsIso_signal[j][iR][1]->Draw("pl same");
			c_effi_frac[j][iR]->SaveAs(Form("pdf_sumIso_combined_effi_frac/pt%d%s_%s.pdf",(int)ptThr,dirName.Data(),c_effi_frac[j][iR]->GetName()));
		}//eta
		c_effi_frac_etaTogether[iR]->SaveAs(Form("pdf_sumIso_combined_effi_frac/pt%d%s_%s.pdf",(int)ptThr,dirName.Data(),c_effi_frac_etaTogether[iR]->GetName()));
	}//radius
}// main

