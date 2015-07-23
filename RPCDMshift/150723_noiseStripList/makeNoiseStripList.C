#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>

void makeNoiseStripList(int run=251562){
	TString fname =Form("DQM_V0001_R000%d__SingleMuon__Run2015B-PromptReco-v1__DQMIO.root",run);
	TFile* f = new TFile( fname.Data() );
	cout << "Input File : " << fname.Data() << endl;
	TH1D* Occ_Barrel=new TH1D("Occ_Barrel","Occupancy of Barrel; Occupancy; Entries", 800,0,8000);
	TH2D* Barrel[5][12]; //[wheel][sector] 
	for(int iwheel=0; iwheel<5; iwheel++){
		for(int isec=0; isec<12; isec++){
			Barrel[iwheel][isec]=(TH2D*)f->Get(Form("DQMData/Run %d/RPC/Run summary/AllHits/Barrel/Wheel_%d/SummaryBySectors/Occupancy_Wheel_%d_Sector_%d",run,iwheel-2,iwheel-2,isec+1));
			//cout << "Occupancy_Wheel_"<<iwheel-2<<"_Sector_"<<isec+1<<endl;
			for(int istrip=0; istrip<Barrel[iwheel][isec]->GetNbinsX();istrip++){
				for(int istage=0; istage<Barrel[iwheel][isec]->GetNbinsY();istage++){
					double occVal = Barrel[iwheel][isec]->GetBinContent(istrip,istage);
					Occ_Barrel->Fill(occVal);
					TString label = Barrel[iwheel][isec]->GetYaxis()->GetBinLabel(istage);	
					if(occVal>1000){
						cout << "W"<<iwheel-2<<"_S"<<isec+1<<" "<<label<<" stripNo."<<istrip<<" (occupancy Value = "<<occVal<<") "<< endl; 	
					}
				}//istage
			}//istrip
		}//isec
	}//iwheel
	
	TH1D* Occ_Endcap=new TH1D("Occ_Endcap","Occupancy of Endcap; Occupancy; Entries", 800,0,8000);
	TH2D* EndcapPlus[4][2];
	TH2D* EndcapMinus[4][2];
	for(int idisk=0; idisk<4; idisk++){
		for(int iring=0; iring<2; iring++){
			TString chMin, chMax;
			if(iring==0) { chMin="CH01"; chMax="CH18";}
			else if(iring==1) { chMin="CH19"; chMax="CH36";}
			EndcapPlus[idisk][iring]=(TH2D*)f->Get(Form("DQMData/Run %d/RPC/Run summary/AllHits/Endcap+/Disk_%d/SummaryByRings/Occupancy_Disk_%d_Ring_%d_%s-%s",run,idisk+1,idisk+1,iring+2,chMin.Data(),chMax.Data()));
			EndcapMinus[idisk][iring]=(TH2D*)f->Get(Form("DQMData/Run %d/RPC/Run summary/AllHits/Endcap-/Disk_%d/SummaryByRings/Occupancy_Disk_%d_Ring_%d_%s-%s",run,-1*(idisk+1),-1*(idisk+1),iring+2,chMin.Data(),chMax.Data()));
			for(int istrip=0; istrip<EndcapPlus[idisk][iring]->GetNbinsX();istrip++){
				for(int istage=0; istage<EndcapPlus[idisk][iring]->GetNbinsY();istage++){
					double occValPlus = EndcapPlus[idisk][iring]->GetBinContent(istrip,istage);
					double occValMinus= EndcapMinus[idisk][iring]->GetBinContent(istrip,istage);
					Occ_Endcap->Fill(occValPlus);
					Occ_Endcap->Fill(occValMinus);
					TString labelPlus = EndcapPlus[idisk][iring]->GetYaxis()->GetBinLabel(istage);	
					TString labelMinus = EndcapMinus[idisk][iring]->GetYaxis()->GetBinLabel(istage);	
					TString roll="RollA";
					int realStripNo=0;
					if(istrip<=32) { roll="RollA"; realStripNo=istrip;}
					else if(istrip>33 && istrip<=32*2) { roll="RollB"; realStripNo=istrip-32;}
					else if(istrip>32*2+1 && istrip<=32*3) { roll="RollC"; realStripNo=istrip-32*2;}

					if(occValPlus>1000)
						cout << "RE+"<<idisk+1<<" "<<labelPlus<<" "<<roll<<" stripNo."<<realStripNo<<" (occupancy Value = "<<occValPlus<<") "<< endl;	
					if(occValMinus>1000)
						cout << "RE"<<-1*(idisk+1)<<" "<<labelMinus<<" "<<roll<<" stripNo."<<realStripNo<<" (occupancy Value = "<<occValMinus<<") "<< endl;
				}//istage
			}//istrip
		}//iring
	}//idisk
	
	TCanvas* c1=new TCanvas();
	Occ_Barrel->Draw();
	c1->SaveAs("Occupancy_of_Barrel.pdf");
	TCanvas* c2=new TCanvas();
	Occ_Endcap->Draw();
	c2->SaveAs("Occupancy_of_Endcap.pdf");
}

