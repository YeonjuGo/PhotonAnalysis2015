/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */

#include "gedPhotonUtility.h" 

const int rocNbins = 100;
const double centBin[nCentCut] = {0,20,60,100,200};
void drawRoc(TGraph*& g1, TH1* h1, TH1* h2)
{
	Int_t n = rocNbins;
	Double_t sigEff[n], bkgEff[n];
	for(int ix=0;ix<n;ix++){
		sigEff[ix] = (h1->Integral(1,ix+1)) / h1->Integral();
		bkgEff[ix] = (h2->Integral(1,ix+1)) / h2->Integral();
		//	cout << "bin " << ix+1 << " >> sigEff : " << sigEff[ix] << ", bkgEffi : " << bkgEff[ix] << endl;
	}
	g1 = new TGraph(rocNbins,sigEff,bkgEff);
	g1->GetXaxis()->SetTitle("N_{sig}[after cut]/N_{sig}[tot]");
	g1->GetYaxis()->SetTitle("N_{bkg}[after cut]/N_{bkg}[tot]");
        g1->SetTitle("ROC curve");
}

void gedPhoton_roc(float ptThr=30, condition cond_ = noC )
{
	gStyle->SetOptStat(0);
	gStyle->SetHistLineWidth(2);

	TString dirName = getCondDirName(cond_);
	TString outSuffix = getCondSuffix(cond_); 

	const char* fName1 = Form("histFiles/gedPhotonHist_AllQCDPhoton30_ggHiNtuplizer_ptThr%d%s.root",(int)ptThr,outSuffix.Data());
	const char* fName2 = Form("histFiles/gedPhotonHist_EmEnrichedDijet30_ggHiNtuplizer_ptThr%d%s.root",(int)ptThr,outSuffix.Data());
	const char* fName3 = Form("histFiles/gedPhotonHist_AllQCDPhoton30_ggHiNtuplizerGED_ptThr%d%s.root",(int)ptThr,outSuffix.Data());
	const char* fName4 = Form("histFiles/gedPhotonHist_EmEnrichedDijet30_ggHiNtuplizerGED_ptThr%d%s.root",(int)ptThr,outSuffix.Data());

	const int NFILE = 4;
	TFile* fin[NFILE];
	fin[0] = new TFile(fName1, "READ");
	fin[1] = new TFile(fName2, "READ");
	fin[2] = new TFile(fName3, "READ");
	fin[3] = new TFile(fName4, "READ");
	std::cout << "input HiForest 1 : " << fin[0]->GetName() << std::endl;
	std::cout << "input HiForest 2 : " << fin[1]->GetName() << std::endl;
	std::cout << "input HiForest 3 : " << fin[2]->GetName() << std::endl;
	std::cout << "input HiForest 4 : " << fin[3]->GetName() << std::endl;

	TH1::SetDefaultSumw2();

	TH1D* h_pho_sumIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_sum_pfIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_sum_pfVsIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_sum_pfBKGIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_ecalIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_pfpVsIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_hcalIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_pfnVsIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_trackerIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2
	TH1D* h_pfcVsIso[nEtaCut][nCentCut][nMomIdCut][NFILE][nRadius];// 2 : file 1 & 2

	for(int j=1; j<nEtaCut; ++j){
		for(int k=1; k<nCentCut; ++k){
			for(int ifile=0; ifile<NFILE; ifile++){
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

					h_ecalIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("pho_ecalClusterIso%s",nameR.Data()));
					h_ecalIso[j][k][0][ifile][iR] -> SetName( Form("pho_ecalClusterIso%s_ifile%d",nameR.Data(),ifile) );
					h_hcalIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("pho_hcalRechitIso%s",nameR.Data()));
					h_hcalIso[j][k][0][ifile][iR] -> SetName( Form("pho_hcalRechitIso%s_ifile%d",nameR.Data(),ifile) );
					h_trackerIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("pho_trackIsoPtCut20%s",nameR.Data()));
					h_trackerIso[j][k][0][ifile][iR] -> SetName( Form("pho_trackIsoPtCut20%s_ifile%d",nameR.Data(),ifile) );
					h_pfpVsIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("pfpVsIso%s",nameR.Data()));
					h_pfpVsIso[j][k][0][ifile][iR] -> SetName( Form("pfpVsIso%s_ifile%d",nameR.Data(),ifile) ); 
					h_pfnVsIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("pfnVsIso%s",nameR.Data()));
					h_pfnVsIso[j][k][0][ifile][iR] -> SetName( Form("pfnVsIso%s_ifile%d",nameR.Data(),ifile) ); 
					h_pfcVsIso[j][k][0][ifile][iR] =  (TH1D*) fin[ifile]-> Get(Form("pfcVsIso%s",nameR.Data()));
					h_pfcVsIso[j][k][0][ifile][iR] -> SetName( Form("pfcVsIso%s_ifile%d",nameR.Data(),ifile) ); 
				}
			}
		}
	}

	TString etaName[4];
	etaName[1] = "|#eta| < 1.4442";
	etaName[2] = "1.566 < |#eta| < 2.0";
	etaName[3] = "|#eta| > 2.0";

	//===================================================
	//===================================================
	// sumIso distribution
#if 0
	TH1D* htemp= new TH1D("htemp",";sumIso;",350/5,-100,250);
	TCanvas* c_sumIso_etaCombine_OLD[nCentCut][nRadius];
	TCanvas* c_sumIso_etaCombine_GED[nCentCut][nRadius];
	for(int iR=2; iR<nRadius; iR++){
		for(int k=1; k<nCentCut; ++k){
			TString name = Form("etaCombine_R%d_cent%d_momId0",iR,k);
			c_sumIso_etaCombine_OLD[k][iR] = new TCanvas(Form("c_sumIso_%s_OLD",name.Data()),"", 1000,300);
			c_sumIso_etaCombine_OLD[k][iR]->Divide(3,1);	
			c_sumIso_etaCombine_GED[k][iR] = new TCanvas(Form("c_sumIso_%s_GED",name.Data()),"", 1000,300);
			c_sumIso_etaCombine_GED[k][iR]->Divide(3,1);	

			for(int j=1; j<nEtaCut; ++j){
				TLegend* leg = new TLegend(0.55,0.65,0.8,0.85);
				leg->SetBorderSize(0);
				leg->SetFillStyle(0);
				TString radiusName = Form("R 0.%d",iR);
				TString centName = Form("%d < hiBin < %d", (int)centBin[k-1], (int)centBin[k]);
				
				c_sumIso_etaCombine_OLD[k][iR]->cd(j);
				htemp->Draw();
				htemp->SetAxisRange(0.00001,0.03,"Y");
				gPad->SetLogy();
				hLineStyle(h_pho_sumIso[j][k][0][0][iR],1,6,1,0,2);
				hLineStyle(h_pho_sumIso[j][k][0][1][iR],1,6,2,0,6);
				h_pho_sumIso[j][k][0][0][iR]->Draw("hist same");
				h_pho_sumIso[j][k][0][1][iR]->Draw("hist same");
				hLineStyle(h_sum_pfVsIso[j][k][0][0][iR],1,8,1,0,3);
				hLineStyle(h_sum_pfVsIso[j][k][0][1][iR],1,8,2,0,8);
				h_sum_pfVsIso[j][k][0][0][iR]->Draw("hist same");
				h_sum_pfVsIso[j][k][0][1][iR]->Draw("hist same");
				leg->AddEntry(h_pho_sumIso[j][k][0][0][iR],"oldIso Sig","l");
				leg->AddEntry(h_pho_sumIso[j][k][0][1][iR],"oldIso Bkg","l");
				leg->AddEntry(h_sum_pfVsIso[j][k][0][0][iR],"pfVsIso Sig","l");
				leg->AddEntry(h_sum_pfVsIso[j][k][0][1][iR],"pfVsIso Bkg","l");
				if(j==3)leg->Draw("same");
				float xp=0.57;
				float yp=0.9;
				float dy=0.06;
				if(j==1){drawText("OLD RECO",xp,yp-3*dy); drawText(radiusName,xp,yp-4*dy); drawText(centName,xp,yp-5*dy);}
				drawText(etaName[j],xp,yp-6*dy);
	
				c_sumIso_etaCombine_GED[k][iR]->cd(j);
				htemp->Draw();
				htemp->SetAxisRange(0.00001,0.03,"Y");
				gPad->SetLogy();
				hLineStyle(h_pho_sumIso[j][k][0][2][iR],1,6,1,0,2);
				hLineStyle(h_pho_sumIso[j][k][0][3][iR],1,6,2,0,6);
				h_pho_sumIso[j][k][0][2][iR]->Draw("hist same");
				h_pho_sumIso[j][k][0][3][iR]->Draw("hist same");
				hLineStyle(h_sum_pfVsIso[j][k][0][2][iR],1,8,1,0,3);
				hLineStyle(h_sum_pfVsIso[j][k][0][3][iR],1,8,2,0,8);
				h_sum_pfVsIso[j][k][0][2][iR]->Draw("hist same");
				h_sum_pfVsIso[j][k][0][3][iR]->Draw("hist same");
				if(j==3)leg->Draw("same");
				if(j==1){drawText("GED RECO",xp,yp-3*dy); drawText(radiusName,xp,yp-4*dy); drawText(centName,xp,yp-5*dy);}
				drawText(etaName[j],xp,yp-6*dy);
			}
			c_sumIso_etaCombine_OLD[k][iR]->SaveAs(Form("pdf_roc%s/pt%d_%s.png",dirName.Data(),(int)ptThr,c_sumIso_etaCombine_OLD[k][iR]->GetName()));
			c_sumIso_etaCombine_GED[k][iR]->SaveAs(Form("pdf_roc%s/pt%d_%s.png",dirName.Data(),(int)ptThr,c_sumIso_etaCombine_GED[k][iR]->GetName()));
		}
	}
#endif
	//===================================================
	//===================================================
	// pfIso background distribution
#if 0
	TH1D* htempBkg= new TH1D("htempBkg",";Iso;",350/5,-100,250);
	TCanvas* c_pfBkg_etaCombine_OLD[nCentCut][nRadius];
	TCanvas* c_pfBkg_etaCombine_GED[nCentCut][nRadius];
	for(int iR=2; iR<nRadius; iR++){
		for(int k=1; k<nCentCut; ++k){
			TString name = Form("etaCombine_R%d_cent%d_momId0",iR,k);
			c_pfBkg_etaCombine_OLD[k][iR] = new TCanvas(Form("c_pfBkg_%s_OLD",name.Data()),"", 1000,300);
			c_pfBkg_etaCombine_OLD[k][iR]->Divide(3,1);	
			c_pfBkg_etaCombine_GED[k][iR] = new TCanvas(Form("c_pfBkg_%s_GED",name.Data()),"", 1000,300);
			c_pfBkg_etaCombine_GED[k][iR]->Divide(3,1);	

			for(int j=1; j<nEtaCut; ++j){
				TLegend* leg = new TLegend(0.55,0.65,0.8,0.85);
				leg->SetBorderSize(0);
				leg->SetFillStyle(0);
				TString radiusName = Form("R 0.%d",iR);
				TString centName = Form("%d < hiBin < %d", (int)centBin[k-1], (int)centBin[k]);

				c_pfBkg_etaCombine_OLD[k][iR]->cd(j);
				htempBkg->Draw();
				htempBkg->SetAxisRange(0.00001,0.03,"Y");
				gPad->SetLogy();		
				hLineStyle(h_sum_pfIso[j][k][0][0][iR],1,1);
				hLineStyle(h_sum_pfVsIso[j][k][0][0][iR],1,2);
				hLineStyle(h_sum_pfBKGIso[j][k][0][0][iR],1,5,1,3001,5);
				h_sum_pfIso[j][k][0][0][iR]->Draw("hist same");
				h_sum_pfVsIso[j][k][0][0][iR]->Draw("hist same");
				h_sum_pfBKGIso[j][k][0][0][iR]->Draw("hist same");
				leg->AddEntry(h_sum_pfIso[j][k][0][0][iR],"pfIso", "l");
				leg->AddEntry(h_sum_pfVsIso[j][k][0][0][iR],"pfVsIso", "l");
				leg->AddEntry(h_sum_pfBKGIso[j][k][0][0][iR],"BkgIso", "l");
				if(j==3)leg->Draw("same");
				float xp=0.5;
				float yp=0.9;
				float dy=0.06;
				if(j==1){drawText("AllQCDPhoton30 OLD RECO",xp,yp-3*dy); drawText(radiusName,xp,yp-4*dy); drawText(centName,xp,yp-5*dy);}
				drawText(etaName[j],xp,yp-6*dy);
	
				c_pfBkg_etaCombine_GED[k][iR]->cd(j);
				htempBkg->Draw();
				gPad->SetLogy();		
				hLineStyle(h_sum_pfIso[j][k][0][2][iR],1,1);
				hLineStyle(h_sum_pfVsIso[j][k][0][2][iR],1,2);
				hLineStyle(h_sum_pfBKGIso[j][k][0][2][iR],1,5,1,3001,5);
				h_sum_pfIso[j][k][0][2][iR]->Draw("hist same");
				h_sum_pfVsIso[j][k][0][2][iR]->Draw("hist same");
				h_sum_pfBKGIso[j][k][0][2][iR]->Draw("hist same");
				if(j==3)leg->Draw("same");
				if(j==1){drawText("AllQCDPhoton30 GED RECO",xp,yp-3*dy); drawText(radiusName,xp,yp-4*dy); drawText(centName,xp,yp-5*dy);}
				drawText(etaName[j],xp,yp-6*dy);
			}
			c_pfBkg_etaCombine_OLD[k][iR]->SaveAs(Form("pdf_roc%s/pt%d_%s.png",dirName.Data(),(int)ptThr,c_pfBkg_etaCombine_OLD[k][iR]->GetName()));
			c_pfBkg_etaCombine_GED[k][iR]->SaveAs(Form("pdf_roc%s/pt%d_%s.png",dirName.Data(),(int)ptThr,c_pfBkg_etaCombine_GED[k][iR]->GetName()));
		}
	}
#endif	

	//===================================================
	//===================================================
	// ROC curve (sumIso)
#if 1
	TGraph* g_roc_sumIso_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfVsIso_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfIso_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_sumIso_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfVsIso_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfIso_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];

	TCanvas* c_roc_etaCombine[nCentCut][nRadius];
	for(int k=1; k<nCentCut; ++k){
		for(int iR=2; iR<nRadius; ++iR){
			TString name = Form("etaCombine_R%d_cent%d_momId0",iR,k);
			c_roc_etaCombine[k][iR] = new TCanvas(Form("c_roc_%s",name.Data()),"", 1000,350);
			c_roc_etaCombine[k][iR]->Divide(3,1);
			TLegend* leg = new TLegend(0.2,0.65,0.7,0.85);
			TString radiusName = Form("R 0.%d",iR);
			for(int j=1; j<nEtaCut; ++j){
				TString centName = Form("%d < hiBin < %d", (int)centBin[k-1], (int)centBin[k]);
				TString nameR = Form("R%d_eta%d_cent%d_momId0",iR,j,k);
				c_roc_etaCombine[k][iR]->cd(j);

				drawRoc(g_roc_sumIso_OLD[j][k][0][iR],h_pho_sumIso[j][k][0][0][iR], h_pho_sumIso[j][k][0][1][iR]);
				//cout << "graph address : " << g_roc_sumIso_OLD[j][k][0][iR] << endl;
				g_roc_sumIso_OLD[j][k][0][iR]->SetName(Form("g_roc_sumIso_OLD_%s",nameR.Data()));
				graphStyle(g_roc_sumIso_OLD[j][k][0][iR], 1, 4);//graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
				g_roc_sumIso_OLD[j][k][0][iR]->Draw("AL");

				drawRoc(g_roc_pfVsIso_OLD[j][k][0][iR],h_sum_pfVsIso[j][k][0][0][iR], h_sum_pfVsIso[j][k][0][1][iR]);
				g_roc_pfVsIso_OLD[j][k][0][iR]->SetName(Form("g_roc_pfVsIso_OLD_%s",nameR.Data()));
				graphStyle(g_roc_pfVsIso_OLD[j][k][0][iR], 2, 4);//graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
				g_roc_pfVsIso_OLD[j][k][0][iR]->Draw("L");

				drawRoc(g_roc_sumIso_GED[j][k][0][iR],h_pho_sumIso[j][k][0][2][iR], h_pho_sumIso[j][k][0][3][iR]);
				g_roc_sumIso_GED[j][k][0][iR]->SetName(Form("g_roc_sumIso_GED_%s",nameR.Data()));
				graphStyle(g_roc_sumIso_GED[j][k][0][iR], 1, 2);//graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
				g_roc_sumIso_GED[j][k][0][iR]->Draw("L");

				drawRoc(g_roc_pfVsIso_GED[j][k][0][iR],h_sum_pfVsIso[j][k][0][2][iR], h_sum_pfVsIso[j][k][0][3][iR]);
				g_roc_pfVsIso_GED[j][k][0][iR]->SetName(Form("g_roc_pfVsIso_GED_%s",nameR.Data()));
				graphStyle(g_roc_pfVsIso_GED[j][k][0][iR], 2, 2);//graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
				g_roc_pfVsIso_GED[j][k][0][iR]->Draw("L");
				float dy=0.06;
				drawText(etaName[j],0.2,0.8-6*dy);
				if(j==1){
					drawText(radiusName,0.2,0.8-4*dy);
					drawText(centName,0.2,0.8-5*dy);
					leg->AddEntry(g_roc_sumIso_OLD[j][k][0][iR], "oldIso OLDreco","l");
					leg->AddEntry(g_roc_pfVsIso_OLD[j][k][0][iR], "pfVsIso OLDreco","l");
					leg->AddEntry(g_roc_sumIso_GED[j][k][0][iR], "oldIso GEDreco","l");
					leg->AddEntry(g_roc_pfVsIso_GED[j][k][0][iR], "pfVsIso GEDreco","l");
					leg->SetBorderSize(0);
					leg->SetFillStyle(0);
					leg->Draw("same");
				}
			}			
			c_roc_etaCombine[k][iR]->SaveAs(Form("pdf_roc%s/pt%d_%s.png",dirName.Data(),(int)ptThr,c_roc_etaCombine[k][iR]->GetName()));
		}
	}
#endif

	//===================================================
	//===================================================
	// ROC curve (Iso ROC seperately)
#if 1
	TGraph* g_roc_oldIso_ecal_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfVsIso_ecal_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_oldIso_ecal_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfVsIso_ecal_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];

	TGraph* g_roc_oldIso_hcal_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfVsIso_hcal_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_oldIso_hcal_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfVsIso_hcal_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];

	TGraph* g_roc_oldIso_tracker_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfVsIso_tracker_OLD[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_oldIso_tracker_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];
	TGraph* g_roc_pfVsIso_tracker_GED[nEtaCut][nCentCut][nMomIdCut][nRadius];

	TCanvas* c_roc_etaCombine_ecal[nCentCut][nRadius];
	TCanvas* c_roc_etaCombine_hcal[nCentCut][nRadius];
	TCanvas* c_roc_etaCombine_tracker[nCentCut][nRadius];
	for(int k=1; k<nCentCut; ++k){
		for(int iR=2; iR<nRadius; ++iR){
			TString name = Form("etaCombine_R%d_cent%d_momId0",iR,k);
			c_roc_etaCombine_ecal[k][iR] = new TCanvas(Form("c_roc_%s_ecal",name.Data()),"", 1000,350);
			c_roc_etaCombine_ecal[k][iR]->Divide(3,1);
			c_roc_etaCombine_hcal[k][iR] = new TCanvas(Form("c_roc_%s_hcal",name.Data()),"", 1000,350);
			c_roc_etaCombine_hcal[k][iR]->Divide(3,1);
			c_roc_etaCombine_tracker[k][iR] = new TCanvas(Form("c_roc_%s_tracker",name.Data()),"", 1000,350);
			c_roc_etaCombine_tracker[k][iR]->Divide(3,1);
			TString radiusName = Form("R 0.%d",iR);
			for(int j=1; j<nEtaCut; ++j){
				TString centName = Form("%d < hiBin < %d", (int)centBin[k-1], (int)centBin[k]);
				TString nameR = Form("R%d_eta%d_cent%d_momId0",iR,j,k);

				c_roc_etaCombine_ecal[k][iR]->cd(j);
				drawRoc(g_roc_oldIso_ecal_OLD[j][k][0][iR],h_ecalIso[j][k][0][0][iR], h_ecalIso[j][k][0][1][iR]);
				g_roc_oldIso_ecal_OLD[j][k][0][iR]->SetName(Form("g_roc_oldIso_ecal_OLD_%s",nameR.Data()));
				graphStyle(g_roc_oldIso_ecal_OLD[j][k][0][iR], 1, 4);//graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
				g_roc_oldIso_ecal_OLD[j][k][0][iR]->Draw("AL");
				drawRoc(g_roc_pfVsIso_ecal_OLD[j][k][0][iR],h_pfpVsIso[j][k][0][0][iR], h_pfpVsIso[j][k][0][1][iR]);
				g_roc_pfVsIso_ecal_OLD[j][k][0][iR]->SetName(Form("g_roc_pfVsIso_ecal_OLD_%s",nameR.Data()));
				graphStyle(g_roc_pfVsIso_ecal_OLD[j][k][0][iR], 2, 4);
				g_roc_pfVsIso_ecal_OLD[j][k][0][iR]->Draw("L");
				drawRoc(g_roc_oldIso_ecal_GED[j][k][0][iR],h_ecalIso[j][k][0][2][iR], h_ecalIso[j][k][0][3][iR]);
				g_roc_oldIso_ecal_GED[j][k][0][iR]->SetName(Form("g_roc_oldIso_ecal_GED_%s",nameR.Data()));
				graphStyle(g_roc_oldIso_ecal_GED[j][k][0][iR], 1, 2);
				g_roc_oldIso_ecal_GED[j][k][0][iR]->Draw("L");
				drawRoc(g_roc_pfVsIso_ecal_GED[j][k][0][iR],h_pfpVsIso[j][k][0][2][iR], h_pfpVsIso[j][k][0][3][iR]);
				g_roc_pfVsIso_ecal_GED[j][k][0][iR]->SetName(Form("g_roc_pfVsIso_ecal_GED_%s",nameR.Data()));
				graphStyle(g_roc_pfVsIso_ecal_GED[j][k][0][iR], 2, 2);
				g_roc_pfVsIso_ecal_GED[j][k][0][iR]->Draw("L");
				float dy=0.06;
				drawText(etaName[j],0.2,0.8-6*dy);
				if(j==1){
					drawText(radiusName,0.2,0.8-4*dy);
					drawText(centName,0.2,0.8-5*dy);
					TLegend* leg = new TLegend(0.2,0.65,0.7,0.85);
					leg->AddEntry(g_roc_oldIso_ecal_OLD[j][k][0][iR], "ecalClusterIso OLD","l");
					leg->AddEntry(g_roc_pfVsIso_ecal_OLD[j][k][0][iR], "pfpVsIso(photon) OLD","l");
					leg->AddEntry(g_roc_oldIso_ecal_GED[j][k][0][iR], "ecalClusterIso GED","l");
					leg->AddEntry(g_roc_pfVsIso_ecal_GED[j][k][0][iR], "pfpVsIso(photon) GED","l");
					leg->SetBorderSize(0);
					leg->SetFillStyle(0);
					leg->Draw("same");
				}

				c_roc_etaCombine_hcal[k][iR]->cd(j);
				drawRoc(g_roc_oldIso_hcal_OLD[j][k][0][iR],h_hcalIso[j][k][0][0][iR], h_hcalIso[j][k][0][1][iR]);
				g_roc_oldIso_hcal_OLD[j][k][0][iR]->SetName(Form("g_roc_oldIso_hcal_OLD_%s",nameR.Data()));
				graphStyle(g_roc_oldIso_hcal_OLD[j][k][0][iR], 1, 4);//graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
				g_roc_oldIso_hcal_OLD[j][k][0][iR]->Draw("AL");
				drawRoc(g_roc_pfVsIso_hcal_OLD[j][k][0][iR],h_pfpVsIso[j][k][0][0][iR], h_pfpVsIso[j][k][0][1][iR]);
				g_roc_pfVsIso_hcal_OLD[j][k][0][iR]->SetName(Form("g_roc_pfVsIso_hcal_OLD_%s",nameR.Data()));
				graphStyle(g_roc_pfVsIso_hcal_OLD[j][k][0][iR], 2, 4);
				g_roc_pfVsIso_hcal_OLD[j][k][0][iR]->Draw("L");
				drawRoc(g_roc_oldIso_hcal_GED[j][k][0][iR],h_hcalIso[j][k][0][2][iR], h_hcalIso[j][k][0][3][iR]);
				g_roc_oldIso_hcal_GED[j][k][0][iR]->SetName(Form("g_roc_oldIso_hcal_GED_%s",nameR.Data()));
				graphStyle(g_roc_oldIso_hcal_GED[j][k][0][iR], 1, 2);
				g_roc_oldIso_hcal_GED[j][k][0][iR]->Draw("L");
				drawRoc(g_roc_pfVsIso_hcal_GED[j][k][0][iR],h_pfpVsIso[j][k][0][2][iR], h_pfpVsIso[j][k][0][3][iR]);
				g_roc_pfVsIso_hcal_GED[j][k][0][iR]->SetName(Form("g_roc_pfVsIso_hcal_GED_%s",nameR.Data()));
				graphStyle(g_roc_pfVsIso_hcal_GED[j][k][0][iR], 2, 2);
				g_roc_pfVsIso_hcal_GED[j][k][0][iR]->Draw("L");
				drawText(etaName[j],0.2,0.8-6*dy);
				if(j==1){
					drawText(radiusName,0.2,0.8-4*dy);
					drawText(centName,0.2,0.8-5*dy);
					TLegend* leg = new TLegend(0.2,0.65,0.7,0.85);
					leg->AddEntry(g_roc_oldIso_hcal_OLD[j][k][0][iR], "hcalRechitIso OLD","l");
					leg->AddEntry(g_roc_pfVsIso_hcal_OLD[j][k][0][iR], "pfnVsIso(neutral) OLD","l");
					leg->AddEntry(g_roc_oldIso_hcal_GED[j][k][0][iR], "hcalRechitIso GED","l");
					leg->AddEntry(g_roc_pfVsIso_hcal_GED[j][k][0][iR], "pfnVsIso(neutral) GED","l");
					leg->SetBorderSize(0);
					leg->SetFillStyle(0);
					leg->Draw("same");
				}

				c_roc_etaCombine_tracker[k][iR]->cd(j);
				drawRoc(g_roc_oldIso_tracker_OLD[j][k][0][iR],h_trackerIso[j][k][0][0][iR], h_trackerIso[j][k][0][1][iR]);
				g_roc_oldIso_tracker_OLD[j][k][0][iR]->SetName(Form("g_roc_oldIso_tracker_OLD_%s",nameR.Data()));
				graphStyle(g_roc_oldIso_tracker_OLD[j][k][0][iR], 1, 4);//graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
				g_roc_oldIso_tracker_OLD[j][k][0][iR]->Draw("AL");
				drawRoc(g_roc_pfVsIso_tracker_OLD[j][k][0][iR],h_pfpVsIso[j][k][0][0][iR], h_pfpVsIso[j][k][0][1][iR]);
				g_roc_pfVsIso_tracker_OLD[j][k][0][iR]->SetName(Form("g_roc_pfVsIso_tracker_OLD_%s",nameR.Data()));
				graphStyle(g_roc_pfVsIso_tracker_OLD[j][k][0][iR], 2, 4);
				g_roc_pfVsIso_tracker_OLD[j][k][0][iR]->Draw("L");
				drawRoc(g_roc_oldIso_tracker_GED[j][k][0][iR],h_trackerIso[j][k][0][2][iR], h_trackerIso[j][k][0][3][iR]);
				g_roc_oldIso_tracker_GED[j][k][0][iR]->SetName(Form("g_roc_oldIso_tracker_GED_%s",nameR.Data()));
				graphStyle(g_roc_oldIso_tracker_GED[j][k][0][iR], 1, 2);
				g_roc_oldIso_tracker_GED[j][k][0][iR]->Draw("L");
				drawRoc(g_roc_pfVsIso_tracker_GED[j][k][0][iR],h_pfpVsIso[j][k][0][2][iR], h_pfpVsIso[j][k][0][3][iR]);
				g_roc_pfVsIso_tracker_GED[j][k][0][iR]->SetName(Form("g_roc_pfVsIso_tracker_GED_%s",nameR.Data()));
				graphStyle(g_roc_pfVsIso_tracker_GED[j][k][0][iR], 2, 2);
				g_roc_pfVsIso_tracker_GED[j][k][0][iR]->Draw("L");
				drawText(etaName[j],0.2,0.8-6*dy);
				if(j==1){
					drawText(radiusName,0.2,0.8-4*dy);
					drawText(centName,0.2,0.8-5*dy);
					TLegend* leg = new TLegend(0.2,0.65,0.7,0.85);
					leg->AddEntry(g_roc_oldIso_tracker_OLD[j][k][0][iR], "trackerIsoPtCut20 OLD","l");
					leg->AddEntry(g_roc_pfVsIso_tracker_OLD[j][k][0][iR], "pfcVsIso(charged) OLD","l");
					leg->AddEntry(g_roc_oldIso_tracker_GED[j][k][0][iR], "trackerIsoPtCut20 GED","l");
					leg->AddEntry(g_roc_pfVsIso_tracker_GED[j][k][0][iR], "pfcVsIso(charged) GED","l");
					leg->SetBorderSize(0);
					leg->SetFillStyle(0);
					leg->Draw("same");
				}



			}			
			c_roc_etaCombine_ecal[k][iR]->SaveAs(Form("pdf_roc%s/pt%d_%s.png",dirName.Data(),(int)ptThr,c_roc_etaCombine_ecal[k][iR]->GetName()));
			c_roc_etaCombine_hcal[k][iR]->SaveAs(Form("pdf_roc%s/pt%d_%s.png",dirName.Data(),(int)ptThr,c_roc_etaCombine_hcal[k][iR]->GetName()));
			c_roc_etaCombine_tracker[k][iR]->SaveAs(Form("pdf_roc%s/pt%d_%s.png",dirName.Data(),(int)ptThr,c_roc_etaCombine_tracker[k][iR]->GetName()));
		}
	}
#endif
}// main

