#!/bin/csh
source /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6/setup.csh; source /afs/cern.ch/work/k/kyolee/public/root_v5.28.00d/bin/thisroot.csh
cp /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/Scripts/pPb_8rap9pt_20150106_rap-1.03--0.43_pT14.0-30.0_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0.csh /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/fit2DData.h /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/fit2DData_all.cpp .
/afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/Fit2DDataAll -f /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_Data_pPb/outRoo_Data_pPb.root 1 -m /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_NPMC_pPb_mcTwoWay/outRoo_NPMC_pPb_mcTwoWay.root /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_PRMC_pPb_mcTwoWay/outRoo_PRMC_pPb_mcTwoWay.root -v sigCB2WNG1 expFunct -d pPb_8rap9pt_20150106 -r etHFm 0 -u 2 -a 1 0 -b 3 1 1 -p 14.0-30.0 -y -1.03--0.43 -t 0.0-100.0 -s 0.000-1.571 -l 1.5-3.0 -x 4 1 /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/outCtErr/fit_ctauErrorRange_pPb.txt -z 0 -j 0 -h 0.0-120.0 -n 0.0-350.0 >& pPb_8rap9pt_20150106_rap-1.03--0.43_pT14.0-30.0_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0.log;
tar zcvf pPb_8rap9pt_20150106_rap-1.03--0.43_pT14.0-30.0_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0.tgz pPb_8rap9pt_20150106_rap-1.03--0.43_pT14.0-30.0_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0* fit2DData.h fit2DData_all.cpp
cp pPb_8rap9pt_20150106_rap-1.03--0.43_pT14.0-30.0_cent0.0-100.0_dPhi0.000-1.571_ntrk0.0-350.0_ET0.0-120.0.tgz /afs/cern.ch/work/k/kyolee/private/CMSSW_5_3_19_Fit/src/pAJpsiFit_mcTwoWay_8rap9pt20150106/Results/pPb_8rap9pt_20150106
