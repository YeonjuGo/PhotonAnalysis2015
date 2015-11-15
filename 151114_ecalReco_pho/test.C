#include "../HiForestAnalysis/hiForest.h"
#include <vector>
void test()
{
   const char* infName = "/afs/cern.ch/work/y/ygo/public/PFphoton/merged_10k_SingleGammaFlatPt10To200_pythia8_FOREST_753p1.root";
   cout << "s" << endl; 
    HiForest h(infName,"",cPbPb);
   h.verbose=0;
//   h.InitTree();
   cout << "s" << endl; 
//   TFile *outf = new TFile("photonTree.root","recreate");
   cout << "s" << endl; 
}

int main(int argc, char** argv)
{
    test();
    return 0;
}
