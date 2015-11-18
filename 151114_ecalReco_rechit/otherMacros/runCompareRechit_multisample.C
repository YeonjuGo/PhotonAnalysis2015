#include "compareRechit_multisample.C"

std::vector<std::string> CollectFiles(TString list) {

  std::vector<std::string> urls;

  ifstream inputFile;
  inputFile.open(list.Data());
  TString line;
  while( inputFile>>line) {
    urls.push_back(line.Data());
  }

  return urls;

}

void runCompareRechit_multisample(TString list = "forestmultifiles/dijet15_global.list", TString list2 = "forestmultifiles/dijet15_multifit.list", const char* cap="DijetMC") {

  list = "forestmultifiles/inputJet80_global.list";
  list2 = "forestmultifiles/jet80_multifit.list";
  std::vector<std::string> urls1 = CollectFiles(list);
  std::vector<std::string> urls2 = CollectFiles(list2);

  //Run analysis
  compareRechit_multisample(urls1, urls2, cap);
 
}

