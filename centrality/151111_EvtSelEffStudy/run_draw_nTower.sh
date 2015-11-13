#void draw_nTower(float etThr=0.0, float eThr=5.0, float norm_i=400, float norm_f=1200, int N=50)

root -l -b -q 'draw_nTower.C+(0.0,1.0)'
#root -l -b -q 'draw_nTower.C+(0.0,1.5)'
root -l -b -q 'draw_nTower.C+(0.0,3.0,400,1000)'
root -l -b -q 'draw_nTower.C+(0.0,5.0,300,800)'
root -l -b -q 'draw_nTower.C+(0.0,7.0,200,600)'
root -l -b -q 'draw_nTower.C+(0.0,9.0,200,600)'

#root -l -b -q 'draw_nTower.C+(0.5,0.0,200,700)'
root -l -b -q 'draw_nTower.C+(1.0,0.0,200,700)'
root -l -b -q 'draw_nTower.C+(2.0,0.0,200,400)'
#root -l -b -q 'draw_nTower.C+(1.5,0.0,200,600)'

#root -l -b -q 'draw_nTower.C+(0.5,1.5,200,600)'
