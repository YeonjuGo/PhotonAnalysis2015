currentDir=`pwd`
#voms-proxy-init --voms cms
echo "srm://cluster142.knu.ac.kr:8443/srm/managerv2?SFN="$1
mkdir -p  ${currentDir}outputDir
var1=$currentDir/outputDir
lcg-ls --offset 0 --count 999 "srm://cluster142.knu.ac.kr:8443/srm/managerv2?SFN="$1 >& temp1.txt
lcg-ls --offset 999 --count 999 "srm://cluster142.knu.ac.kr:8443/srm/managerv2?SFN="$1 >& temp1.txt
cat temp1.txt | awk -v p1=$var1  '{ print "srmcp -debug=true -srm_protocol_version=2 srm://cluster142.knu.ac.kr:8443/srm/managerv2?SFN="$1" file:///" p1$1 }' |     bash

#rm temp1.txt
