
# void leadingPho_draw_JEC(const char *calgo="ak3PF", bool savePlots=1, bool fitting=1)

# use gaus fitting 

root -l -b -q 'leadingPho_draw_JEC.C++("ak3PF","ak3PF_gaus",1)'
root -l -b -q 'leadingPho_draw_JEC.C++("ak4PF","ak4PF_gaus",1)'
root -l -b -q 'leadingPho_draw_JEC.C++("ak5PF","ak5PF_gaus",1)'

root -l -b -q 'leadingPho_draw_JEC.C++("akPu3PF","akPu3PF_gaus",1)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu4PF","akPu4PF_gaus",1)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu5PF","akPu5PF_gaus",1)'

:<<"end"
root -l -b -q 'leadingPho_draw_JEC.C++("ak3Calo",1)'
root -l -b -q 'leadingPho_draw_JEC.C++("ak4Calo",1)'
root -l -b -q 'leadingPho_draw_JEC.C++("ak5Calo",1)'

root -l -b -q 'leadingPho_draw_JEC.C++("akPu3Calo",1)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu4Calo",1)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu5Calo",1)'
end

# use mean value
root -l -b -q 'leadingPho_draw_JEC.C++("ak3PF","ak3PF_mean",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("ak4PF","ak4PF_mean",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("ak5PF","ak5PF_mean",1,0)'

root -l -b -q 'leadingPho_draw_JEC.C++("akPu3PF","akPu3PF_mean",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu4PF","akPu4PF_mean",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu5PF","akPu5PF_mean",1,0)'

:<< 'sss'
root -l -b -q 'leadingPho_draw_JEC.C++("ak3Calo",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("ak4Calo",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("ak5Calo",1,0)'

root -l -b -q 'leadingPho_draw_JEC.C++("akPu3Calo",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu4Calo",1,0)'
root -l -b -q 'leadingPho_draw_JEC.C++("akPu5Calo",1,0)'
sss
