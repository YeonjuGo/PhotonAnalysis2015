#include "TFile.h"
#include "TTree.h"
#include "yjUtility.h"
#include <iostream>
using namespace std;
const int colhere[] = {1,2,4,9,6,46,kOrange};

TH1D* GetHisthiBin(string run="695");
const char* dir = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIMinimumBias2/Merged/";

void temp_compare_runs_1D_hiBin(){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();

    const int Nrun = 4;
    string run[] = {"620","656","694","816"};
    //string run[] = {"695","548","620","656","694","816"};
    //string run[] = {"695","620","656","694"};
    const char* cap = "";


    TH1D* h[Nrun];
    TLegend* l1 = new TLegend(0.4,0.6,0.95,0.9);
    legStyle(l1);
    TCanvas* c1 = new TCanvas("c1","", 500, 500);
    for(int i=0; i<Nrun; ++i){
        h[i]=GetHisthiBin(run[i]);
        SetHistColor(h[i],col[i]);
        if(i==0) h[i]->DrawCopy("hist");
        else h[i]->DrawCopy("hist same");
        l1->AddEntry(h[i], Form("Run 262%s",run[i].data()));

    }
    l1->Draw();
    c1->SaveAs(Form("pdf/compareBtwRuns_run262%s_hiBin%s_offlineCentBin.pdf",run[0].data(),cap));

    TCanvas* c2 = new TCanvas("c2","", 500, 500);
    for(int i=1; i<Nrun; ++i){
        h[i]->Divide(h[0]);
        h[i]->SetMarkerStyle(20);
        h[i]->SetMarkerSize(0.85);
        h[i]->GetYaxis()->SetRangeUser(0.7,1.3);
        h[i]->SetTitle(Form(";hiBin;Run XXX / Run 262695"));
        if(i==1) h[i]->DrawCopy("pe");
        else h[i]->DrawCopy("pe same");
    }
    jumSun(0,1,200,1);
    c2->SaveAs(Form("pdf/compareBtwRuns_ratio_run262%s_hiBin%s_offlineCentBin.pdf",run[0].data(),cap));


}

TH1D* GetHisthiBin(string run){
    TFile* f1;
    if(run=="811") f1 = TFile::Open(Form("/afs/cern.ch/work/y/ygo/public/PFphoton/HIForestMinimumBias_run262811_file1.root"));
    else if(run=="816") f1 = TFile::Open(Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIMinimumBias_run262816.root"));
    else if(run=="548") f1 = TFile::Open(Form("root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HIForestMinbiasUPC_run262548.root"));
    else f1 = TFile::Open(Form("%sHIForestExpress_run262%s.root",dir,run.data()));
    //if(f1->IsZombie()) { cout << run.data() << " doesn't exist!! " << endl; continue;}


    float bounds[200] = {
        0, 0, 4.89949, 7.99844, 8.92695, 9.66829, 10.4113, 11.0529, 11.7822, 12.4533, 13.1472, 13.9668, 14.7675, 15.5483, 16.4161, 17.3047, 18.1283, 19.0679, 20.051, 21.0405, 22.2351, 23.2812, 24.5173, 25.897, 27.1474, 28.459, 29.9013, 31.3157, 32.7494, 34.3873, 36.0111, 37.7366, 39.4742, 41.2605, 43.0352, 45.1332, 47.3762, 49.6964, 51.879, 54.4101, 56.8655, 59.4357, 61.9746, 64.6677, 67.6169, 70.5101, 73.4362, 76.9011, 80.2528, 83.2244, 86.7075, 90.251, 94.0596, 98.2499, 102.294, 106.763, 111.263, 116.182, 120.897, 125.363, 130.158, 135.085, 140.409, 146.08, 151.909, 157.415, 163.238, 169.32, 175.464, 182.866, 189.667, 197.613, 204.473, 212.603, 220.045, 227.635, 235.659, 244.188, 253.487, 261.769, 269.849, 278.664, 288.425, 297.848, 306.746, 316.694, 327.771, 337.893, 348.05, 359.054, 370.201, 381.906, 394.094, 406.909, 419.204, 432.827, 446.383, 460.286, 474.495, 489.801, 503.344, 517.537, 531.617, 545.604, 562.185, 577.949, 594.962, 612.415, 630.409, 647.859, 665.238, 683.566, 701.988, 721.231, 740.548, 759.612, 778.229, 797.684, 819.26, 844.711, 865.76, 888.009, 910.645, 931.656, 952.614, 974.873, 1001.99, 1026.11, 1050.47, 1074.6, 1096.31, 1123.83, 1148.94, 1174.01, 1199.23, 1228.52, 1258.23, 1286.83, 1317.29, 1347.72, 1376.62, 1405.23, 1440.22, 1472.31, 1505, 1539.24, 1572.41, 1607.52, 1638.3, 1675.36, 1714.35, 1749.18, 1787.1, 1826.05, 1866.69, 1905.91, 1944.2, 1985.24, 2025.87, 2069.84, 2116.89, 2163.97, 2207.1, 2253.35, 2297.07, 2343.68, 2386.1, 2436.58, 2485.76, 2539.44, 2590.07, 2641.06, 2692.61, 2740.54, 2797.53, 2853.44, 2906.16, 2970.53, 3026.77, 3083.32, 3140.47, 3207.33, 3268.97, 3335.02, 3411.45, 3479.12, 3553.01, 3626.72, 3694.18, 3765.61, 3840.09, 3916.15, 4000.97, 4081.75, 4173.33, 4264.96, 4367.08, 4463.28, 4572.44, 4702.18
    };

    TTree* tref = (TTree*)f1->Get("hiEvtAnalyzer/HiTree");
    float hf; int hiBin;
    tref->SetBranchAddress("hiHF", &hf);
    tref->SetBranchAddress("hiBin", &hiBin);

    TTree *tskim = (TTree*)f1->Get("skimanalysis/HltTree");
    int phfCoincFilter3, pprimaryVertexFilter;
    tskim->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3);
    tskim->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);

    TTree *thlt = (TTree*)f1->Get("hltanalysis/HltTree");
    int hlt_mbhfand;
    thlt->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1", &hlt_mbhfand);

    TFile* outf = new TFile(Form("skimfiles/compare_centralitybins%s.root",run.data()),"recreate");
    int newbin;
    TTree* t = new TTree("anaCentrality","analysis level centrality");
    t->Branch("newBin",&newbin,"newBin/I");
    t->Branch("oldBin",&hiBin,"newBin/I");

    TH1D* h = new TH1D(Form("h%s",run.data()), Form(";hiBin;EventFraction"), 200,0,200);
    int N = tref->GetEntries();
    int nSelEvt = 0;
    for(int i = 0; i < N; ++i){
        tref->GetEntry(i);
        thlt->GetEntry(i);
        tskim->GetEntry(i);
        if(i % 10000 == 0) cout<<"processing event : "<<i<<endl;
        if(phfCoincFilter3 != 1) continue;
        if(pprimaryVertexFilter != 1) continue;
        if(hlt_mbhfand != 1) continue;

        nSelEvt++;
        newbin = 199;
        for(int b = 0; b < 200; ++b){
            if(hf >= bounds[199-b]){
                newbin = b;
                h->Fill(newbin);
                break;	    
            }
        }
        t->Fill();
    }
    h->Scale(1./nSelEvt);
    t->Write(); 
    outf->Write(); 
    return h;
}


