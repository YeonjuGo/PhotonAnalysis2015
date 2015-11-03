/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */

#include "gedPhotonUtility.h" 
void gedPhoton_effi_figures(TString treePath="ggHiNtuplizer", float ptThr=30, condition cond_ = noC )
{
	gStyle->SetOptStat(0);
	gStyle->SetHistLineWidth(2);
	
	TString recoName = "";
	if(treePath=="ggHiNtuplizer") recoName = "OLD";
	else if(treePath=="ggHiNtuplizerGED") recoName = "GED";

        TString dirName = getCondDirName(cond_);
        TString outSuffix = getCondSuffix(cond_); 
    
	const char* fName1 = Form("histFiles/gedPhotonHist_AllQCDPhoton30_%s_ptThr%d%s.root",treePath.Data(),(int)ptThr,outSuffix.Data());
	const char* fName2 = Form("histFiles/gedPhotonHist_EmEnrichedDijet30_%s_ptThr%d%s.root",treePath.Data(),(int)ptThr,outSuffix.Data());
	const char* fNameOut = Form("histFiles/gedPhotonHist_effi_%s_ptThr%d%s.root",treePath.Data(),(int)ptThr,outSuffix.Data());

	TFile* fin[2];
	fin[0] = new TFile(fName1, "READ");
	fin[1] = new TFile(fName2, "READ");
	std::cout << "input HiForest 1 : " << fin[0]->GetName() << std::endl;
	std::cout << "input HiForest 2 : " << fin[1]->GetName() << std::endl;
	
	TFile* fout = new TFile(fNameOut, "RECREATE");

	TH1::SetDefaultSumw2();

	TH1D* h_pho_sumIso[nEtaCut][nCentCut][nMomIdCut][2][nRadius];// 2 : file 1 & 2
	TH1D* h_sum_pfIso[nEtaCut][nCentCut][nMomIdCut][2][nRadius];// 2 : file 1 & 2
	TH1D* h_sum_pfVsIso[nEtaCut][nCentCut][nMomIdCut][2][nRadius];// 2 : file 1 & 2
	TH1D* h_sum_pfBKGIso[nEtaCut][nCentCut][nMomIdCut][2][nRadius];// 2 : file 1 & 2
	for(int j=1; j<nEtaCut; ++j){
		for(int k=1; k<nCentCut; ++k){
			for(int ifile=0; ifile<2; ifile++){
				for(int iR=2; iR<nRadius; iR++){
					TString nameR = Form("R%d_eta%d_cent%d_momId0",iR,j,k);
					h_pho_sumIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("pho_sumIso%s",nameR.Data()));
					h_pho_sumIso[j][k][0][ifile][iR] -> SetName(Form("pho_sumIso%s_ifile%d",nameR.Data(),ifile)); 
					h_sum_pfIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("sum_pfIso%s",nameR.Data()));
					h_sum_pfIso[j][k][0][ifile][iR] -> SetName(Form("sum_pfIso%s_ifile%d",nameR.Data(),ifile)); 
					h_sum_pfVsIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("sum_pfVsIso%s",nameR.Data()));
					h_sum_pfVsIso[j][k][0][ifile][iR] -> SetName(Form("sum_pfVsIso%s_ifile%d",nameR.Data(),ifile)); 
					h_sum_pfBKGIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("sum_pfBKGIso%s",nameR.Data()));
					h_sum_pfBKGIso[j][k][0][ifile][iR] -> SetName(Form("sum_pfBKGIso%s_ifile%d",nameR.Data(),ifile)); 
				}
			}
		}
	}

	double centBin[nCentCut] = {0,20,60,100,200};
	TCanvas* c_sumIso[nEtaCut][nCentCut][nRadius];
	TCanvas* c_pfIso[nEtaCut][nCentCut][nRadius];
	TCanvas* c_pfVsIso[nEtaCut][nCentCut][nRadius];
	TCanvas* c_pfBKGIso[nEtaCut][nCentCut][nRadius];
	TCanvas* c_sumIso_etaCombine[nCentCut][nRadius];
	TCanvas* c_pfIso_etaCombine[nCentCut][nRadius];
	TCanvas* c_pfVsIso_etaCombine[nCentCut][nRadius];
	TCanvas* c_pfBKGIso_etaCombine[nCentCut][nRadius];
	double d_pho_sumIso_effi[nEtaCut][nCentCut][nRadius];
	double d_pho_sumIso_frac[nEtaCut][nCentCut][nRadius];
	double d_pho_sumIso_effi_err[nEtaCut][nCentCut][nRadius];
	double d_pho_sumIso_frac_err[nEtaCut][nCentCut][nRadius];
	double d_pho_sumIso_cross[nEtaCut][nCentCut][nRadius];

	double d_sum_pfIso_effi[nEtaCut][nCentCut][nRadius];
	double d_sum_pfIso_frac[nEtaCut][nCentCut][nRadius];
	double d_sum_pfIso_effi_err[nEtaCut][nCentCut][nRadius];
	double d_sum_pfIso_frac_err[nEtaCut][nCentCut][nRadius];
	double d_sum_pfIso_cross[nEtaCut][nCentCut][nRadius];

	double d_sum_pfVsIso_effi[nEtaCut][nCentCut][nRadius];
	double d_sum_pfVsIso_frac[nEtaCut][nCentCut][nRadius];
	double d_sum_pfVsIso_effi_err[nEtaCut][nCentCut][nRadius];
	double d_sum_pfVsIso_frac_err[nEtaCut][nCentCut][nRadius];
	double d_sum_pfVsIso_cross[nEtaCut][nCentCut][nRadius];

	for(int j=1; j<nEtaCut; ++j){
		for(int k=1; k<nCentCut; ++k){
			for(int iR=2; iR<nRadius; ++iR){
				d_pho_sumIso_effi[j][k][iR]=0.0;
				d_pho_sumIso_frac[j][k][iR]=0.0;
				d_pho_sumIso_effi_err[j][k][iR]=0.0;
				d_pho_sumIso_frac_err[j][k][iR]=0.0;
				d_pho_sumIso_cross[j][k][iR]=0.0;
				d_sum_pfIso_effi[j][k][iR]=0.0;
				d_sum_pfIso_frac[j][k][iR]=0.0;
				d_sum_pfIso_effi_err[j][k][iR]=0.0;
				d_sum_pfIso_frac_err[j][k][iR]=0.0;
				d_sum_pfIso_cross[j][k][iR]=0.0;
				d_sum_pfVsIso_effi[j][k][iR]=0.0;
				d_sum_pfVsIso_frac[j][k][iR]=0.0;
				d_sum_pfVsIso_effi_err[j][k][iR]=0.0;
				d_sum_pfVsIso_frac_err[j][k][iR]=0.0;
				d_sum_pfVsIso_cross[j][k][iR]=0.0;
			}
		}
	}

	for(int k=1; k<nCentCut; ++k){
		for(int iR=2; iR<nRadius; ++iR){
			TString name = Form("etaCombine_R%d_cent%d_momId0",iR,k);
			c_sumIso_etaCombine[k][iR] = new TCanvas(Form("c_sumIso%s",name.Data()),"", 1000,300);
			c_pfIso_etaCombine[k][iR] = new TCanvas(Form("c_sum_pfIso%s",name.Data()),"", 1000,300);
			c_pfVsIso_etaCombine[k][iR] = new TCanvas(Form("c_sum_pfVsIso%s",name.Data()),"", 1000,300);
			c_pfBKGIso_etaCombine[k][iR] = new TCanvas(Form("c_sum_pfBKGIso%s",name.Data()),"", 1000,300);
			c_sumIso_etaCombine[k][iR]->Divide(3,1);
			c_pfIso_etaCombine[k][iR]->Divide(3,1);
			c_pfVsIso_etaCombine[k][iR]->Divide(3,1);
			c_pfBKGIso_etaCombine[k][iR]->Divide(3,1);

			for(int j=1; j<nEtaCut; ++j){
				TString centName = Form("%d < hiBin < %d", (int)centBin[k-1], (int)centBin[k]);
				TString nameR = Form("R%d_eta%d_cent%d_momId0",iR,j,k);
				c_sumIso[j][k][iR] = new TCanvas(Form("c_sumIso%s",nameR.Data()),"", 800,600);
				c_pfIso[j][k][iR] = new TCanvas(Form("c_sum_pfIso%s",nameR.Data()),"", 800,600);
				c_pfVsIso[j][k][iR] = new TCanvas(Form("c_sum_pfVsIso%s",nameR.Data()),"", 800,600);
				c_pfBKGIso[j][k][iR] = new TCanvas(Form("c_sum_pfBKGIso%s",nameR.Data()),"", 800,600);

				c_sumIso[j][k][iR]->cd();
				h_pho_sumIso[j][k][0][0][iR]->Draw("hist");
				h_pho_sumIso[j][k][0][1][iR]->SetLineColor(2);
				h_pho_sumIso[j][k][0][1][iR]->Draw("hist same");

				c_pfIso[j][k][iR]->cd();
				h_sum_pfIso[j][k][0][0][iR]->Draw("hist");
				h_sum_pfIso[j][k][0][1][iR]->SetLineColor(2);
				h_sum_pfIso[j][k][0][1][iR]->Draw("hist same");

				c_pfVsIso[j][k][iR]->cd();
				h_sum_pfVsIso[j][k][0][0][iR]->Draw("hist");
				h_sum_pfVsIso[j][k][0][1][iR]->SetLineColor(2);
				h_sum_pfVsIso[j][k][0][1][iR]->Draw("hist same");

				c_pfBKGIso[j][k][iR]->cd();
				h_sum_pfBKGIso[j][k][0][0][iR]->Draw("hist");
				h_sum_pfBKGIso[j][k][0][1][iR]->SetLineColor(2);
				h_sum_pfBKGIso[j][k][0][1][iR]->Draw("hist same");

				c_sumIso_etaCombine[k][iR]->cd(j);
				h_pho_sumIso[j][k][0][0][iR]->Draw("hist");
				h_pho_sumIso[j][k][0][1][iR]->SetLineColor(2);
				h_pho_sumIso[j][k][0][1][iR]->Draw("hist same");

				c_pfIso_etaCombine[k][iR]->cd(j);
				h_sum_pfIso[j][k][0][0][iR]->Draw("hist");
				h_sum_pfIso[j][k][0][1][iR]->SetLineColor(2);
				h_sum_pfIso[j][k][0][1][iR]->Draw("hist same");

				c_pfVsIso_etaCombine[k][iR]->cd(j);
				h_sum_pfVsIso[j][k][0][0][iR]->Draw("hist");
				h_sum_pfVsIso[j][k][0][1][iR]->SetLineColor(2);
				h_sum_pfVsIso[j][k][0][1][iR]->Draw("hist same");

				c_pfBKGIso_etaCombine[k][iR]->cd(j);
				h_sum_pfBKGIso[j][k][0][0][iR]->Draw("hist");
				h_sum_pfBKGIso[j][k][0][1][iR]->SetLineColor(2);
				h_sum_pfBKGIso[j][k][0][1][iR]->Draw("hist same");


			
				//==========================
				//==========================
				// frac & effi calculation for old calo 
				// error calculation
				// poisson error for fraction , binomial error for efficiency.
				c_sumIso_etaCombine[k][iR]->cd(j);
				double tmpCross = 0;
				tmpCross = findCross(h_pho_sumIso[j][k][0][0][iR], h_pho_sumIso[j][k][0][1][iR], d_pho_sumIso_frac[j][k][iR], d_pho_sumIso_effi[j][k][iR], d_pho_sumIso_frac_err[j][k][iR], d_pho_sumIso_effi_err[j][k][iR]);
				d_pho_sumIso_cross[j][k][iR] = tmpCross;
				float dy=0.06;
				drawText(Form("cross value : %.1f", d_pho_sumIso_cross[j][k][iR]) ,0.5,0.8);
				drawText(Form("purity : %.3f", d_pho_sumIso_frac[j][k][iR]) ,0.5,0.8-dy);
				drawText(Form("efficiency : %.3f", d_pho_sumIso_effi[j][k][iR]) ,0.5,0.8-2*dy);
				drawText(Form("eff*pur : %.3f", d_pho_sumIso_effi[j][k][iR]*d_pho_sumIso_frac[j][k][iR]) ,0.5,0.8-3*dy);
				drawText(recoName ,0.5,0.8-4*dy);
				drawText(centName ,0.6,0.8-4*dy);
				gPad->SetLogy();
				jumSun(tmpCross,0,tmpCross,h_pho_sumIso[j][k][0][0][iR]->GetMaximum());

				c_pfIso_etaCombine[k][iR]->cd(j);
				tmpCross = 0;
				tmpCross = findCross(h_sum_pfIso[j][k][0][0][iR], h_sum_pfIso[j][k][0][1][iR], d_sum_pfIso_frac[j][k][iR], d_sum_pfIso_effi[j][k][iR], d_sum_pfIso_frac_err[j][k][iR], d_sum_pfIso_effi_err[j][k][iR]);
				d_sum_pfIso_cross[j][k][iR] = tmpCross;
				drawText(Form("cross value : %.1f", d_sum_pfIso_cross[j][k][iR]) ,0.5,0.8);
				drawText(Form("purity : %.3f", d_sum_pfIso_frac[j][k][iR]) ,0.5,0.8-dy);
				drawText(Form("efficiency : %.3f", d_sum_pfIso_effi[j][k][iR]) ,0.5,0.8-2*dy);
				drawText(Form("eff*pur : %.3f", d_sum_pfIso_effi[j][k][iR]*d_sum_pfIso_frac[j][k][iR]) ,0.5,0.8-3*dy);
				drawText(recoName ,0.5,0.8-4*dy);
				drawText(centName ,0.6,0.8-4*dy);
				gPad->SetLogy();
				jumSun(tmpCross,0,tmpCross,h_sum_pfIso[j][k][0][0][iR]->GetMaximum());
			
				c_pfVsIso_etaCombine[k][iR]->cd(j);
				tmpCross = 0;
				tmpCross = findCross(h_sum_pfVsIso[j][k][0][0][iR], h_sum_pfVsIso[j][k][0][1][iR], d_sum_pfVsIso_frac[j][k][iR], d_sum_pfVsIso_effi[j][k][iR], d_sum_pfVsIso_frac_err[j][k][iR], d_sum_pfVsIso_effi_err[j][k][iR]);
				d_sum_pfVsIso_cross[j][k][iR] = tmpCross;
				drawText(Form("cross value : %.1f", d_sum_pfVsIso_cross[j][k][iR]) ,0.5,0.8);
				drawText(Form("purity : %.3f", d_sum_pfVsIso_frac[j][k][iR]) ,0.5,0.8-dy);
				drawText(Form("efficiency : %.3f", d_sum_pfVsIso_effi[j][k][iR]) ,0.5,0.8-2*dy);
				drawText(Form("eff*pur : %.3f", d_sum_pfVsIso_effi[j][k][iR]*d_sum_pfVsIso_frac[j][k][iR]) ,0.5,0.8-3*dy);
				drawText(recoName ,0.5,0.8-4*dy);
				drawText(centName ,0.6,0.8-4*dy);
				gPad->SetLogy();
				jumSun(tmpCross,0,tmpCross,h_sum_pfVsIso[j][k][0][0][iR]->GetMaximum());
		
				c_sumIso[j][k][iR]->cd();
				tmpCross = 0;
				tmpCross = findCross(h_pho_sumIso[j][k][0][0][iR], h_pho_sumIso[j][k][0][1][iR], d_pho_sumIso_frac[j][k][iR], d_pho_sumIso_effi[j][k][iR], d_pho_sumIso_frac_err[j][k][iR], d_pho_sumIso_effi_err[j][k][iR]);
				d_pho_sumIso_cross[j][k][iR] = tmpCross;
				dy=0.06;
				drawText(Form("cross value : %.1f", d_pho_sumIso_cross[j][k][iR]) ,0.5,0.8);
				drawText(Form("purity : %.3f", d_pho_sumIso_frac[j][k][iR]) ,0.5,0.8-dy);
				drawText(Form("efficiency : %.3f", d_pho_sumIso_effi[j][k][iR]) ,0.5,0.8-2*dy);
				drawText(Form("eff*pur : %.3f", d_pho_sumIso_effi[j][k][iR]*d_pho_sumIso_frac[j][k][iR]) ,0.5,0.8-3*dy);
				drawText(recoName ,0.5,0.8-4*dy);
				drawText(centName ,0.6,0.8-4*dy);
				gPad->SetLogy();
				jumSun(tmpCross,0,tmpCross,h_pho_sumIso[j][k][0][0][iR]->GetMaximum());
				c_sumIso[j][k][iR]->SaveAs(Form("pdf_sum_detIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_sumIso[j][k][iR]->GetName()));

				c_pfIso[j][k][iR]->cd();
				tmpCross = 0;
				tmpCross = findCross(h_sum_pfIso[j][k][0][0][iR], h_sum_pfIso[j][k][0][1][iR], d_sum_pfIso_frac[j][k][iR], d_sum_pfIso_effi[j][k][iR], d_sum_pfIso_frac_err[j][k][iR], d_sum_pfIso_effi_err[j][k][iR]);
				d_sum_pfIso_cross[j][k][iR] = tmpCross;
				drawText(Form("cross value : %.1f", d_sum_pfIso_cross[j][k][iR]) ,0.5,0.8);
				drawText(Form("purity : %.3f", d_sum_pfIso_frac[j][k][iR]) ,0.5,0.8-dy);
				drawText(Form("efficiency : %.3f", d_sum_pfIso_effi[j][k][iR]) ,0.5,0.8-2*dy);
				drawText(Form("eff*pur : %.3f", d_sum_pfIso_effi[j][k][iR]*d_sum_pfIso_frac[j][k][iR]) ,0.5,0.8-3*dy);
				drawText(recoName ,0.5,0.8-4*dy);
				drawText(centName ,0.6,0.8-4*dy);
				gPad->SetLogy();
				jumSun(tmpCross,0,tmpCross,h_sum_pfIso[j][k][0][0][iR]->GetMaximum());
				c_pfIso[j][k][iR]->SaveAs(Form("pdf_sum_pfIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_pfIso[j][k][iR]->GetName()));

				c_pfVsIso[j][k][iR]->cd();
				tmpCross = 0;
				tmpCross = findCross(h_sum_pfVsIso[j][k][0][0][iR], h_sum_pfVsIso[j][k][0][1][iR], d_sum_pfVsIso_frac[j][k][iR], d_sum_pfVsIso_effi[j][k][iR], d_sum_pfVsIso_frac_err[j][k][iR], d_sum_pfVsIso_effi_err[j][k][iR]);
				d_sum_pfVsIso_cross[j][k][iR] = tmpCross;
				drawText(Form("cross value : %.1f", d_sum_pfVsIso_cross[j][k][iR]) ,0.5,0.8);
				drawText(Form("purity : %.3f", d_sum_pfVsIso_frac[j][k][iR]) ,0.5,0.8-dy);
				drawText(Form("efficiency : %.3f", d_sum_pfVsIso_effi[j][k][iR]) ,0.5,0.8-2*dy);
				drawText(Form("eff*pur : %.3f", d_sum_pfVsIso_effi[j][k][iR]*d_sum_pfVsIso_frac[j][k][iR]) ,0.5,0.8-3*dy);
				drawText(recoName ,0.5,0.8-4*dy);
				drawText(centName ,0.6,0.8-4*dy);
				gPad->SetLogy();
				jumSun(tmpCross,0,tmpCross,h_sum_pfVsIso[j][k][0][0][iR]->GetMaximum());
				c_pfVsIso[j][k][iR]->SaveAs(Form("pdf_sum_pfVsIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_pfVsIso[j][k][iR]->GetName()));
			}			

			c_pfVsIso_etaCombine[k][iR]->SaveAs(Form("pdf_sum_pfVsIso_%s%s/pt%d_etaCombine_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_pfVsIso_etaCombine[k][iR]->GetName()));
			c_pfIso_etaCombine[k][iR]->SaveAs(Form("pdf_sum_pfIso_%s%s/pt%d_etaCombine_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_pfIso_etaCombine[k][iR]->GetName()));
			c_sumIso_etaCombine[k][iR]->SaveAs(Form("pdf_sum_detIso_%s%s/pt%d_etaCombine_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_sumIso_etaCombine[k][iR]->GetName()));
		}
	}

	//=================================================================================
	//=================================================================================
	// centrality dependent efficiency & S/B ratio histogram. 

	TH1D* h_detIso_effi[nEtaCut][nRadius]; //as a fnc. of hiBin
	TH1D* h_detIso_frac[nEtaCut][nRadius]; //as a fnc. of hiBin
	TH1D* h_detIso_signal[nEtaCut][nRadius]; //as a fnc. of hiBin
	TH1D* h_pfIso_effi[nEtaCut][nRadius]; //as a fnc. of hiBin
	TH1D* h_pfIso_frac[nEtaCut][nRadius]; //as a fnc. of hiBin
	TH1D* h_pfIso_signal[nEtaCut][nRadius]; //as a fnc. of hiBin
	TH1D* h_pfVsIso_effi[nEtaCut][nRadius]; //as a fnc. of hiBin
	TH1D* h_pfVsIso_frac[nEtaCut][nRadius]; //as a fnc. of hiBin
	TH1D* h_pfVsIso_signal[nEtaCut][nRadius]; //as a fnc. of hiBin
	TCanvas* c_detIso_effi[nEtaCut][nRadius];
	TCanvas* c_detIso_frac[nEtaCut][nRadius];
	TCanvas* c_detIso_signal[nEtaCut][nRadius];
	TCanvas* c_pfIso_effi[nEtaCut][nRadius];
	TCanvas* c_pfIso_frac[nEtaCut][nRadius];	
	TCanvas* c_pfIso_signal[nEtaCut][nRadius];	
	TCanvas* c_pfVsIso_effi[nEtaCut][nRadius];
	TCanvas* c_pfVsIso_frac[nEtaCut][nRadius];
	TCanvas* c_pfVsIso_signal[nEtaCut][nRadius];

	fout -> cd();
	for(int j=1; j<nEtaCut; ++j){
		for(int iR=2; iR<nRadius; ++iR){
			TString name = Form("R%d_eta%d_momId0",iR,j);
			h_detIso_effi[j][iR] = new TH1D(Form("h_detIso_effi_%s",name.Data()), "sumIso efficiency ;hiBin;Effficiency", nCentCut-1, centBin); 
			h_detIso_frac[j][iR] = new TH1D(Form("h_detIso_frac_%s",name.Data()), "sumIso purity ;hiBin;Purity", nCentCut-1, centBin);	
			h_detIso_signal[j][iR] = new TH1D(Form("h_detIso_signal_%s",name.Data()), "sumIso effi*pur ;hiBin;Efficiency*Purity", nCentCut-1, centBin);	
			h_pfIso_effi[j][iR] = new TH1D(Form("h_pfIso_effi_%s",name.Data()), "sumIso efficiency ;hiBin;Effficiency", nCentCut-1, centBin); 
			h_pfIso_frac[j][iR] = new TH1D(Form("h_pfIso_frac_%s",name.Data()), "sumIso purity ;hiBin;Purity", nCentCut-1, centBin);	
			h_pfIso_signal[j][iR] = new TH1D(Form("h_pfIso_signal_%s",name.Data()), "sumIso effi*pur ;hiBin;Efficiency*Purity", nCentCut-1, centBin);	
			h_pfVsIso_effi[j][iR] = new TH1D(Form("h_pfVsIso_effi_%s",name.Data()), "sumIso efficiency ;hiBin;Effficiency", nCentCut-1, centBin); 
			h_pfVsIso_frac[j][iR] = new TH1D(Form("h_pfVsIso_frac_%s",name.Data()), "sumIso purity ;hiBin;Purity", nCentCut-1, centBin);	
			h_pfVsIso_signal[j][iR] = new TH1D(Form("h_pfVsIso_signal_%s",name.Data()), "sumIso effi*pur ;hiBin;Efficiency*Purity", nCentCut-1, centBin);	
			for(int k=1; k<nCentCut; ++k){
				h_detIso_effi[j][iR]->SetBinContent(k, d_pho_sumIso_effi[j][k][iR]);
				h_detIso_frac[j][iR]->SetBinContent(k, d_pho_sumIso_frac[j][k][iR]);
				h_detIso_signal[j][iR]->SetBinContent(k, d_pho_sumIso_effi[j][k][iR] * d_pho_sumIso_frac[j][k][iR]);
			//	h_detIso_effi[j][iR]->SetBinError(k, d_pho_sumIso_effi_err[j][k][iR]);
			//	h_detIso_frac[j][iR]->SetBinError(k, d_pho_sumIso_frac_err[j][k][iR]);
		//		cout << "d_pho_sumIso_effi[j][k][iR] : " << d_pho_sumIso_effi[j][k][iR] << endl;
		//		cout << "d_pho_sumIso_frac[j][k][iR] : " << d_pho_sumIso_frac[j][k][iR] << endl;
		//		cout << "d_pho_sumIso_effi_err[j][k][iR] : " << d_pho_sumIso_effi_err[j][k][iR] << endl;
		//		cout << "d_pho_sumIso_frac_err[j][k][iR] : " << d_pho_sumIso_frac_err[j][k][iR] << endl;
				h_pfIso_effi[j][iR]->SetBinContent(k, d_sum_pfIso_effi[j][k][iR]);
				h_pfIso_frac[j][iR]->SetBinContent(k, d_sum_pfIso_frac[j][k][iR]);
				h_pfIso_signal[j][iR]->SetBinContent(k, d_sum_pfIso_effi[j][k][iR] * d_sum_pfIso_frac[j][k][iR]);
		//		h_pfIso_effi[j][iR]->SetBinError(k, d_sum_pfIso_effi_err[j][k][iR]);
		//		h_pfIso_frac[j][iR]->SetBinError(k, d_sum_pfIso_frac_err[j][k][iR]);
	
				h_pfVsIso_effi[j][iR]->SetBinContent(k, d_sum_pfVsIso_effi[j][k][iR]);
				h_pfVsIso_frac[j][iR]->SetBinContent(k, d_sum_pfVsIso_frac[j][k][iR]);
				h_pfVsIso_signal[j][iR]->SetBinContent(k, d_sum_pfVsIso_effi[j][k][iR] * d_sum_pfVsIso_frac[j][k][iR]);
		//		h_pfVsIso_effi[j][iR]->SetBinError(k, d_sum_pfVsIso_effi_err[j][k][iR]);
		//		h_pfVsIso_frac[j][iR]->SetBinError(k, d_sum_pfVsIso_frac_err[j][k][iR]);
			}//cent

			c_detIso_effi[j][iR] = new TCanvas(Form("c_detIso_effi_%s",name.Data()),"", 800,600);
			c_detIso_effi[j][iR]->cd();
			h_detIso_effi[j][iR]->Draw("hist e");
			c_detIso_frac[j][iR] = new TCanvas(Form("c_detIso_frac_%s",name.Data()),"", 800,600);
			c_detIso_frac[j][iR]->cd();
			h_detIso_frac[j][iR]->Draw("hist e");
			c_pfIso_effi[j][iR] = new TCanvas(Form("c_pfIso_effi_%s",name.Data()),"", 800,600);
			c_pfIso_effi[j][iR]->cd();
			h_pfIso_effi[j][iR]->Draw("hist e");
			c_pfIso_frac[j][iR] = new TCanvas(Form("c_pfIso_frac_%s",name.Data()),"", 800,600);
			c_pfIso_frac[j][iR]->cd();
			h_pfIso_frac[j][iR]->Draw("hist e");
			c_pfVsIso_effi[j][iR] = new TCanvas(Form("c_pfVsIso_effi_%s",name.Data()),"", 800,600);
			c_pfVsIso_effi[j][iR]->cd();
			h_pfVsIso_effi[j][iR]->Draw("hist e");
			c_pfVsIso_frac[j][iR] = new TCanvas(Form("c_pfVsIso_frac_%s",name.Data()),"", 800,600);
			c_pfVsIso_frac[j][iR]->cd();
			h_pfVsIso_frac[j][iR]->Draw("hist e");
			
			///////////////////////////////////
			//save all the histograms and canvas
			h_detIso_effi[j][iR]->Write();	
			h_detIso_frac[j][iR]->Write();	
			h_detIso_signal[j][iR]->Write();	
			h_pfIso_effi[j][iR]->Write();	
			h_pfIso_frac[j][iR]->Write();	
			h_pfIso_signal[j][iR]->Write();	
			h_pfVsIso_effi[j][iR]->Write();	
			h_pfVsIso_frac[j][iR]->Write();	
			h_pfVsIso_signal[j][iR]->Write();	
			c_detIso_effi[j][iR]->SaveAs(Form("pdf_sum_detIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_detIso_effi[j][iR]->GetName()));
			c_detIso_frac[j][iR]->SaveAs(Form("pdf_sum_detIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_detIso_frac[j][iR]->GetName()));
			c_pfIso_effi[j][iR]->SaveAs(Form("pdf_sum_pfIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_detIso_effi[j][iR]->GetName()));
			c_pfIso_frac[j][iR]->SaveAs(Form("pdf_sum_pfIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_detIso_frac[j][iR]->GetName()));
			c_pfVsIso_effi[j][iR]->SaveAs(Form("pdf_sum_pfVsIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_detIso_effi[j][iR]->GetName()));
			c_pfVsIso_frac[j][iR]->SaveAs(Form("pdf_sum_pfVsIso_%s%s/pt%d_%s.pdf",recoName.Data(),dirName.Data(),(int)ptThr,c_detIso_frac[j][iR]->GetName()));
		}//radius
	}//eta	
	fout -> Close();
}// main

