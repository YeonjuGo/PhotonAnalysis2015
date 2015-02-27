#!/bin/csh
source /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6/setup.csh; source /afs/cern.ch/work/k/kyolee/public/root_v5.28.00d/bin/thisroot.csh
cp /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/Scripts/Pbp_8rap9pt_20150106_rap-0.47-0.43_pT7.5-8.5_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0.csh /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/fit2DData.h /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/fit2DData_all.cpp .
/afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/Fit2DDataAll -q 1 /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_Data_Pbp_v2/outRoo_Data_Pbp_v2.root -f /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_Data_Pbp_v1/outRoo_Data_Pbp_v1.root 1 -m /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_NPMC_Pbp_mcTwoWay/outRoo_NPMC_Pbp_mcTwoWay.root /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_PRMC_Pbp_mcTwoWay/outRoo_PRMC_Pbp_mcTwoWay.root -v sigCB2WNG1 expFunct -d Pbp_8rap9pt_20150106 -r etHFm 0 -u 2 -a 1 0 -b 2 1 1 -p 7.5-8.5 -y -0.47-0.43 -t 0.0-100.0 -s 0.000-1.571 -l 1.5-3.0 -x 4 1 /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/outCtErr/fit_ctauErrorRange_Pbp.txt -z 0 -j 0 -h 0.0-120.0 -n 0.0-350.0 >& Pbp_8rap9pt_20150106_rap-0.47-0.43_pT7.5-8.5_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0.log;
tar zcvf Pbp_8rap9pt_20150106_rap-0.47-0.43_pT7.5-8.5_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0.tgz Pbp_8rap9pt_20150106_rap-0.47-0.43_pT7.5-8.5_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0* fit2DData.h fit2DData_all.cpp
cp Pbp_8rap9pt_20150106_rap-0.47-0.43_pT7.5-8.5_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0.tgz /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/Results/Pbp_8rap9pt_20150106
