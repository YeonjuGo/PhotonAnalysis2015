#!/bin/bash
srmls srm://cluster142.knu.ac.kr:8443/srm/managerv2?SFN=/pnfs/knu.ac.kr/data/cms/store/user/ygo/photon_pPb_JEC_forest_15 >& cplist.txt
awk '{print $2}' cplist.txt > tmplist.txt
exec < tmplist.txt
set i=1
while read line
do
	if [ $i -gt 1 ]; then
		root -l -q -b 'forest2yskim_jetSkim_forestV3.C+("/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton120.root", "", 170, 2517,35, "yskimmedFiles/yskim_pA_AllQCDPhoton120170.root", 5,0)'
		srmcp srm://cluster142.knu.ac.kr:8443/srm/managerv2?SFN=$line file:////u/user/goyeonju/PRODUCTION/CMSSW_5_3_20/src/pPb_JEC/merging-forest/pPbAllQCD15
	fi
	i=$(expr $i + 1);
done
