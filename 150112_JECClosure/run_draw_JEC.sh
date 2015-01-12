
#void leadingPho_draw_JEC(
#                const char* infile="merged_old_pPb.root", 
#                        TString calgo="ak3PF", bool savePlots=1, bool gausfitting=1)

# use gaus fitting 

root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak3PF",1,1)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak4PF",1,1)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak5PF",1,1)'

root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu3PF",1,1)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu4PF",1,1)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu5PF",1,1)'

# use mean value
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak3PF",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak4PF",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak5PF",1,0)'

root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu3PF",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu4PF",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu5PF",1,0)'

# Calo algorithms
:<<"end"
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak3Calo",1,1)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak4Calo",1,1)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak5Calo",1,1)'

root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu3Calo",1,1)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu4Calo",1,1)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu5Calo",1,1)'
end
:<< 'sss'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak3Calo",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak4Calo",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","ak5Calo",1,0)'

root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu3Calo",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("merged_old_pPb.root","akPu4Calo",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu5Calo",1,0)'
sss
