/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#include "../gedPhotonUtility.h" 
static const long MAXTREESIZE = 10000000000;

void pfVsSumIsoCutTuning(TString var="pfVs", float ptThr=20)
{

    const char* fName[2];
    fName[0] = "/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_AllQCDPhoton30_PhotonFilter20GeV_eta24_TuneZ2_PbPb_5020GeV.root";
    fName[1] = "/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_EmEnrichedDijet30_ParticleFilter20GeV_eta24_TuneZ2_5020GeV.root";
    TFile* f[2];//[signal, background]
    TTree* t_pho_old[2];//[signal, background]
    TTree* t_pho_ged[2];//[signal, background]
    TTree* t_hiEvt[2];//[signal, background]
    for(int i=0; i<2; i++){
        f[i]= new TFile(fName[i],"READ");
        std::cout << "input HiForest : " << f[i]->GetName() << std::endl;
        t_hiEvt[i] = (TTree*)f[i]->Get("hiEvtAnalyzer/HiTree");
        t_pho_old[i] = (TTree*)f[i]->Get("ggHiNtuplizer/EventTree");
        t_pho_ged[i] = (TTree*)f[i]->Get("ggHiNtuplizerGED/EventTree");
    }

    TH1D* h1D_iso_old[2];
    TH1D* h1D_iso_ged[2];

    for(int i=0; i<2; i++){
        h1D_iso_old[i] = new TH1D(Form("h1D_pfIso_old_isBkg%d",i),";pfVsSumIso4;",300,-50,150);
        h1D_iso_ged[i] = (TH1D*) h1D_iso_old[i]->Clone(Form("h1D_pfIso_ged_isBkg%d",i));
    }

    for(int i=0; i<2; i++){
        t_pho_old[i]->Draw(Form("pfcVsIso4+pfnVsIso4+pfpVsIso4>>%s",h1D_iso_old[i]->GetName()));
        t_pho_ged[i]->Draw(Form("pfcVsIso4+pfnVsIso4+pfpVsIso4>>%s",h1D_iso_ged[i]->GetName()));
        h1D_iso_old[i]->SetLineColor(i+1);
        h1D_iso_ged[i]->SetLineColor(i+1);
    }

    TCanvas* c1= new TCanvas("c1","",1000,300);
    c1->Divide(2,1);
    c1->cd(1);
    h1D_iso_old[0]->Draw("hist");
    h1D_iso_old[1]->Draw("hist same");
    c1->cd(2);
    h1D_iso_ged[0]->Draw("hist");
    h1D_iso_ged[1]->Draw("hist same");

}
