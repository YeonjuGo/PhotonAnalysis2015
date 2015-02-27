#!/bin/csh
set TOP="$PWD"

cd /afs/cern.ch/work/y/ygo/private/PRODUCTION/CMSSW_5_3_20/src/
eval `scramv1 runtime -csh`
cd $TOP
S2_AccDistMaker_batch 1 1 0

