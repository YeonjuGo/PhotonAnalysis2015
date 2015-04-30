#void compare_pf_eff(const int pfid=2, const int pdgid=22, const int ptThr=10, const int centLow = 0, const int centUp=10, const float etaLow = 0, const float etaUp = 1.44)

#electron barrel
root -l -b -q 'compare_pf_eff.C++(2,11,10,0,20,0,1.44)'
root -l -b -q 'compare_pf_eff.C++(2,11,10,20,60,0,1.44)'
root -l -b -q 'compare_pf_eff.C++(2,11,10,60,200,0,1.44)'

#electron endcap
root -l -b -q 'compare_pf_eff.C++(2,11,10,0,20,1.44,5.0)'
root -l -b -q 'compare_pf_eff.C++(2,11,10,20,60,1.44,5.0)'
root -l -b -q 'compare_pf_eff.C++(2,11,10,60,200,1.44,5.0)'

#photon barrel
#root -l -b -q 'compare_pf_eff.C++(4,22,10,0,20,0,1.44)'
#root -l -b -q 'compare_pf_eff.C++(4,22,10,20,60,0,1.44)'
#root -l -b -q 'compare_pf_eff.C++(4,22,10,60,199,0,1.44)'

#photon endcap
#root -l -b -q 'compare_pf_eff.C++(4,22,10,0,20,1.44,5.0)'
#root -l -b -q 'compare_pf_eff.C++(4,22,10,20,60,1.44,5.0)'
#root -l -b -q 'compare_pf_eff.C++(4,22,10,60,199,1.44,5.0)'
