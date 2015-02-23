#void yj_rootAna_acceptance(char *strBinnig = "8rap9pt", bool isPrompt=true, bool isPbp=true)

root -l -b -q 'AccDistMaker.C++("8rap9pt", 0, 0)'
root -l -b -q 'AccDistMaker.C++("8rap9pt", 0, 1)'
root -l -b -q 'AccDistMaker.C++("8rap9pt", 1, 0)'
root -l -b -q 'AccDistMaker.C++("8rap9pt", 1, 1)'

