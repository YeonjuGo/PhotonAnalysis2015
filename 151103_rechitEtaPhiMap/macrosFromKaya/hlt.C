#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include "EventMatchingCMS.h"

#define MAXPHOTONS 500

int hlt() {
    TFile* hltfile = new TFile("/export/d00/scratch/cfmcginn/hltDir/openHLT_20151021_HIGamma30502_740F.root");
    TFile* forestfile = new TFile("/mnt/hadoop/cms/store/user/luck/L1Emulator/Pyquen_Unquenched_AllQCDPhoton30_PhotonFilter20GeV_eta24_TuneZ2_PbPb_5020GeV_HiForest.root");

    TTree* hlttree = (TTree*)(hltfile->Get("hltbitanalysis/HltTree"));
    TTree* foresttree = (TTree*)(forestfile->Get("multiPhotonAnalyzer/photon"));
    TTree* evttree = (TTree*)(forestfile->Get("hiEvtAnalyzer/HiTree"));

    int L1_SingleEG5_BptxAND, L1_SingleEG12_BptxAND, L1_SingleEG15_BptxAND, L1_SingleEG18_BptxAND, L1_SingleEG15_Centrality50_100;
    ULong64_t event;
    int lumi, run;

    hlttree->SetBranchAddress("Event", &event);
    hlttree->SetBranchAddress("LumiBlock", &lumi);
    hlttree->SetBranchAddress("Run", &run);
    hlttree->SetBranchAddress("L1_SingleEG5_BptxAND", &L1_SingleEG5_BptxAND);
    hlttree->SetBranchAddress("L1_SingleEG12_BptxAND", &L1_SingleEG12_BptxAND);
    hlttree->SetBranchAddress("L1_SingleEG15_BptxAND", &L1_SingleEG15_BptxAND);
    hlttree->SetBranchAddress("L1_SingleEG18_BptxAND", &L1_SingleEG18_BptxAND);
    hlttree->SetBranchAddress("L1_SingleEG15_Centrality50_100", &L1_SingleEG15_Centrality50_100);

    Int_t forest_event;
    int forest_lumi, forest_run;
    Int_t nPho;
    Float_t et[MAXPHOTONS];
    Float_t eta[MAXPHOTONS];
    Float_t phi[MAXPHOTONS];
    Float_t cc4[MAXPHOTONS];
    Float_t cr4[MAXPHOTONS];
    Float_t ct4[MAXPHOTONS];
    Float_t r9[MAXPHOTONS];
    Float_t hoverE[MAXPHOTONS];
    Float_t sigmaIetaIeta[MAXPHOTONS];
    Float_t swissCrx[MAXPHOTONS];
    Float_t seedTime[MAXPHOTONS];

    evttree->SetBranchAddress("evt", &forest_event);
    evttree->SetBranchAddress("lumi", &forest_lumi);
    evttree->SetBranchAddress("run", &forest_run);
    foresttree->SetBranchAddress("nPhotons", &nPho);
    foresttree->SetBranchAddress("pt", et);
    foresttree->SetBranchAddress("eta", eta);
    foresttree->SetBranchAddress("phi", phi);
    foresttree->SetBranchAddress("cc4", cc4);
    foresttree->SetBranchAddress("cr4", cr4);
    foresttree->SetBranchAddress("ct4PtCut20", ct4);
    foresttree->SetBranchAddress("r9", r9);
    foresttree->SetBranchAddress("hadronicOverEm", hoverE);
    foresttree->SetBranchAddress("sigmaIetaIeta", sigmaIetaIeta);
    foresttree->SetBranchAddress("swissCrx", swissCrx);
    foresttree->SetBranchAddress("seedTime", seedTime);

    EventMatchingCMS* ematcher = new EventMatchingCMS();

    int matched_entry = -1;
    int hltentries = hlttree->GetEntries();
    for (int i=0; i<hltentries; i++) {
        hlttree->GetEntry(i);
        ematcher->addEvent(event, lumi, run, i);
    }

    TH1D* hpt = new TH1D("hpt", "", 100, 0, 100);
    TH1D* heta = new TH1D("heta", "", 100, -3, 3);
    TH1D* hphi = new TH1D("hphi", "", 100, -4, 4);
    TH1D* hcc4 = new TH1D("hcc4", "", 100, -60, 60);
    TH1D* hcr4 = new TH1D("hcr4", "", 100, -40, 60);
    TH1D* hct4 = new TH1D("hct4", "", 100, -40, 60);
    TH1D* hr9 = new TH1D("hr9", "", 100, 0, 1);
    TH1D* hhoe = new TH1D("hhoe", "", 100, 0, 1);
    TH1D* hsii = new TH1D("hsii", "", 100, 0, 0.1);
    TH2D* hscrx = new TH2D("hscrx", "", 100, -15, 15, 100, -2, 2);

    TH1D* heg5pt = new TH1D("heg5pt", "", 100, 0, 100);
    TH1D* heg12pt = new TH1D("heg12pt", "", 100, 0, 100);
    TH1D* heg15pt = new TH1D("heg15pt", "", 100, 0, 100);
    TH1D* heg18pt = new TH1D("heg18pt", "", 100, 0, 100);
    TH1D* heg15cpt = new TH1D("heg15cpt", "", 100, 0, 100);

    int forestentries = foresttree->GetEntries();
    for (int i=0; i<forestentries; i++) {
        foresttree->GetEntry(i);
        evttree->GetEntry(i);
        matched_entry = ematcher->retrieveEvent(forest_event, forest_lumi, forest_run);
        if (matched_entry == -1)
            continue;

        hlttree->GetEntry(matched_entry);

        hpt->Fill(et[0]);
        heta->Fill(eta[0]);
        hphi->Fill(phi[0]);
        hcc4->Fill(cc4[0]);
        hcr4->Fill(cr4[0]);
        hct4->Fill(ct4[0]);
        hr9->Fill(r9[0]);
        hhoe->Fill(hoverE[0]);
        hsii->Fill(sigmaIetaIeta[0]);
        hscrx->Fill(seedTime[0], swissCrx[0]);

        // if (L1_SingleEG5_BptxAND)
        //     heg5pt->Fill(et[0]);
        // if (L1_SingleEG12_BptxAND)
        //     heg12pt->Fill(et[0]);
        // if (L1_SingleEG15_BptxAND)
        //     heg15pt->Fill(et[0]);
        // if (L1_SingleEG18_BptxAND)
        //     heg18pt->Fill(et[0]);
        // if (L1_SingleEG15_Centrality50_100)
        //     heg15cpt->Fill(et[0]);

        // if (!L1_SingleEG5_BptxAND && et[0] > 15)
        //     printf("L1_SingleEG5_BptxAND  - forestindex: %i - event: %i; pt: %f\n", i, forest_event, et[0]);
        // if (!L1_SingleEG12_BptxAND && et[0] > 22)
        //     printf("L1_SingleEG12_BptxAND - forestindex: %i - event: %i; pt: %f\n", i, forest_event, et[0]);
        // if (!L1_SingleEG15_BptxAND && et[0] > 25)
        //     printf("L1_SingleEG15_BptxAND - forestindex: %i - event: %i; pt: %f\n", i, forest_event, et[0]);
        // if (!L1_SingleEG18_BptxAND && et[0] > 28)
        //     printf("L1_SingleEG18_BptxAND - forestindex: %i - event: %i; pt: %f\n", i, forest_event, et[0]);
        // if (!L1_SingleEG15_Centrality50_100 && et[0] > 100)
        //     printf("L1_SingleEG15_Centrality50_100 - event: %i; pt: %f\n", forest_event, et[0]);

        et[0] = -99;
        eta[0] = -99;
        phi[0] = -99;
        cc4[0] = -99;
        cr4[0] = -99;
        ct4[0] = -99;
        r9[0] = -99;
        hoverE[0] = -99;
        sigmaIetaIeta[0] = -99;
        swissCrx[0] = -99;
        seedTime[0] = -99;

        matched_entry = -1;
        forest_run = 0;
        forest_lumi = 0;
        forest_event = 0;
        run = 0;
        lumi = 0;
        event = 0;
    }

    TCanvas* c1 = new TCanvas("c1", "", 1024, 1024);
    c1->Divide(3, 3);
    c1->cd(1);
    heta->Draw();
    c1->cd(2);
    hphi->Draw();
    c1->cd(3);
    hcc4->Draw();
    c1->cd(4);
    hcr4->Draw();
    c1->cd(5);
    hct4->Draw();
    c1->cd(6);
    hr9->Draw();
    c1->cd(7);
    hhoe->Draw();
    c1->cd(8);
    hsii->Draw();
    c1->cd(9);
    hscrx->Draw("colz");

    TCanvas* c2 = new TCanvas("c2", "", 540, 540);
    hpt->Draw();

    // TCanvas* c1 = new TCanvas("c1", "", 1200, 720);
    // c1->Divide(3, 2);
    // c1->cd(1);
    // heg5pt->Divide(hpt);
    // heg5pt->SetMaximum(1.01);
    // heg5pt->SetMinimum(0.98);
    // heg5pt->Draw();
    // c1->cd(2);
    // heg12pt->Divide(hpt);
    // heg12pt->SetMaximum(1.01);
    // heg12pt->SetMinimum(0.98);
    // heg12pt->Draw();
    // c1->cd(3);
    // heg15pt->Divide(hpt);
    // heg15pt->SetMaximum(1.01);
    // heg15pt->SetMinimum(0.98);
    // heg15pt->Draw();
    // c1->cd(4);
    // heg18pt->Divide(hpt);
    // heg18pt->SetMaximum(1.01);
    // heg18pt->SetMinimum(0.98);
    // heg18pt->Draw();
    // c1->cd(5);
    // heg15cpt->Divide(hpt);
    // heg15cpt->SetMaximum(1.01);
    // heg15cpt->SetMinimum(0.98);
    // heg15cpt->Draw();
    // c1->cd(6);
    // gPad->SetLogy();
    // hpt->Draw();

    return 0;
}
