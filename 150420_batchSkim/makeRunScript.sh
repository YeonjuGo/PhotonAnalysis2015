#!/bin/bash
if [ $# -ne 1 ]; then
#  echo "Usage: $0 [Executable] [Dataset] [Prefix] "
  echo "Usage: $0 [Dataset] "
  exit;
fi

executable=$(pwd)/forest2yskim_jetSkim_forestV3.C
datasets=$1
#prefix=$3

################################################################
########## Script parameter setting
################################################################

storage=$(pwd)/yskimmedFiles/$datasets
if [ ! -d "$(pwd)/yskimmedFiles" ]; then
  mkdir $(pwd)/yskimmedFiles
fi
if [ ! -d "$storage" ]; then
  mkdir $storage
fi

txtrst=$(tput sgr0)
txtred=$(tput setaf 2)  #1 Red, 2 Green, 3 Yellow, 4 Blue, 5 Purple, 6 Cyan, 7 White
txtbld=$(tput bold)     #dim (Half-bright mode), bold (Bold characters)

################################################################
########## Information print out
################################################################

echo "Running macro with ${txtbld}${txtred}$executable${txtrst} on ${txtbld}${txtred}$datasets${txtrst}."
if [ "$storage" != "" ]; then
  echo "yskimmedFiles will be placed on ${txtbld}${txtred}$storage${txtrst}."
fi

script=$datasets"_runForest2yskim_jetSkim_forestV3"

echo "Processing: "$script
printf "#!/bin/bash\n" > $(pwd)/$script.sh
printf "cmsenv5320\n" >> $(pwd)/$script.sh

srmls srm://cluster142.knu.ac.kr:8443/srm/managerv2?SFN=/pnfs/knu.ac.kr/data/cms/store/user/ygo/2014-photon-forests/$datasets >& tmp_list_$datasets.txt

awk '/HiForest/' tmp_list_$datasets.txt > list_$datasets.txt
rm tmp_list_$datasets.txt

while read A B C ETC
do
    stringZ=${B}
    tmpfn=${stringZ:92}
    command="root -l -b -q 'forest2yskim_jetSkim_forestV3.C++(\"dcap://cluster142.knu.ac.kr/${B}\", \"\", 40, \"$storage/yskimmed_$tmpfn\",4,0)'"
    echo $command >> $(pwd)/$script.sh
done < list_$datasets.txt
