#!/bin/csh
set TOP="$PWD"

source /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6/setup.csh;
cd /afs/cern.ch/work/y/ygo/private/PRODUCTION/CMSSW_5_3_20/src/
eval `scramv1 runtime -csh`
cd $TOP
/afs/cern.ch/work/y/ygo/private/PhotonAnalysis2015/150227_JPsiAccSys2ndBatch/S2_AccDistMaker_batch 1 1 2 >& S2_AccDistMaker_batch_1_1_2.log
mv S2_AccDistMaker_batch_1_1_2.log LOG

