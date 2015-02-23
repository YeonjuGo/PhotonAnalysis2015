# void FitAndMakeRatioToy_pt(int nToy=1, bool isPrompt=true, bool isPbp=true)

root -l -b -q 'FitAndMakeRatioToy_pt.C++(100, 0, 0)' 
root -l -b -q 'FitAndMakeRatioToy_pt.C++(100, 0, 1)'
#root -l -b -q 'FitAndMakeRatioToy_pt.C++(100, 1, 0)'
#root -l -b -q 'FitAndMakeRatioToy_pt.C++(100, 1, 1)'

#root -l -b -q 'FitAndMakeRatioToy_pt.C++(100, 0, 0)' > FitAndMakeRatioToy_pt00.log &
#root -l -b -q 'FitAndMakeRatioToy_pt.C++(100, 0, 1)' > FitAndMakeRatioToy_pt01.log &
#root -l -b -q 'FitAndMakeRatioToy_pt.C++(100, 1, 0)' > FitAndMakeRatioToy_pt10.log &
#root -l -b -q 'FitAndMakeRatioToy_pt.C++(100, 1, 1)' > FitAndMakeRatioToy_pt11.log &

