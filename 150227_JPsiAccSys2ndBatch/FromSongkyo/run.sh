#!/bin/bash


#################################################
############ for Pbp (1st run)
#################################################

### for v1 and v2 merging
./runBatch_Pbp.sh Fit2DDataAll /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_Data_Pbp_v1/outRoo_Data_Pbp_v1.root /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_Data_Pbp_v2/outRoo_Data_Pbp_v2.root Pbp_8rap9pt_20150106 >& log_Pbp_8rap9pt_20150106 &

#################################################
############ for pPb (2nd run) 
#################################################
./runBatch_pPb.sh Fit2DDataAll /afs/cern.ch/work/k/kyolee/private/pAJpsi_rooDataSet_mcTwoWay_zVtxCutWeight/outRoo_Data_pPb/outRoo_Data_pPb.root pPb_8rap9pt_20150106 >& log_pPb_8rap9pt_20150106 &

