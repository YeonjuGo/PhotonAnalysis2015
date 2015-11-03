#void draw_gedPhoton_effi(float ptThr=30, condition cond_ = noC )

root -l -b -q 'draw_gedPhoton_effi.C++(30, noC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(30, hoeC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(30, sigmaC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(30, hoeAndSigmaC)'

root -l -b -q 'draw_gedPhoton_effi.C+(40, noC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(40, hoeC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(40, sigmaC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(40, hoeAndSigmaC)'

root -l -b -q 'draw_gedPhoton_effi.C+(60, noC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(60, hoeC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(60, sigmaC)'
#root -l -b -q 'draw_gedPhoton_effi.C+(60, hoeAndSigmaC)'
