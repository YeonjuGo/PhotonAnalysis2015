#!/bin/sh

progName="gedPhotonMacros_RecHit";

g++ $progName.C $(root-config --cflags --libs) -Wall -O2 -o $progName.exe || exit 1

inputFile="";
##inputFile="/mnt/hadoop/cms/store/user/luck/HIMinBiasUPC/2011_MB_750_hiForest/0.root"; # NOTE : this is a big file. GetEntries() = 1662337
#inputFile="/mnt/hadoop/cms/store/user/tatar/HIHighPt/HiForest_HIHighPt_photon30_HIRun2011-v1.root"; # GetEntries() = 166993

##outputFile="gedPhotonMacros_RecHit2_2011_MB_750_hiForest.root";
##outputFile="gedPhotonMacros_PfCand_2011_MB_750_hiForest.root";
###outputFile="gedPhotonMacros_RecHit2_HiForest_HIHighPt_photon30_HIRun2011-v1.root";
#outputFile="gedPhotonMacros_PfCand_HiForest_HIHighPt_photon30_HIRun2011-v1.root";
outputFile="";

./$progName.exe $inputFile $outputFile || exit 1
