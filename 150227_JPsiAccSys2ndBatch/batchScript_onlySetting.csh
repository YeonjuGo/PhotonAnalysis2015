#!/bin/csh
set TOP="$PWD"

source /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6/setup.csh;
cd /afs/cern.ch/work/y/ygo/private/PRODUCTION/CMSSW_5_3_20/src/
eval `scramv1 runtime -csh`
cd $TOP
#S2_AccDistMaker_batch 1 1 0 >& S2_AccDistMaker_batch_1_1_0.log

