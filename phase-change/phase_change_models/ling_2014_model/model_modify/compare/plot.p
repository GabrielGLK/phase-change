
reset
set size square
set title 'Stefan problem interface position'
set xlabel 'Time t(s)'
set ylabel 'delta(mm)'
set key right bottom
p [0:10][0:0.0027]'in-pos-1' u ($1+0.29):2 w p ls 1  t 'Ling-model-1',\
'in-pos-2' u ($1+0.29):2 w p ls 2 t 'Ling-model-2',\
'out-type1' u 1:3 w l lw 2 lc rgb 'red' t 'Analytical solution'
pause mouse key

reset
set size square
set title 'Stefan problem liquid velocity'
set xlabel 'Time t(s)'
set ylabel 'u_l(m/s)'
set key right top
p [0:10][0:0.002]'out-type1' every 20 u ($1+0.29):6 w p ls 1  t 'Ling-model-1',\
'out-type2' every 20 u ($1+0.29):6 w p ls 2  t 'Ling-model-2',\
'out-type1' u 1:5 w l lw 2 lc rgb 'red' t 'Analytical solution'
pause mouse key

// temperature
reset
set size square
set title 'Temperature at t=3s'
set xlabel 'X(mm)'
set ylabel 'Temperature(K)'
set key right top
p [0:0.01]'temperature-1-3.0' u 1:2 w p pt 'o'  t 'Ling-model-1',\
'temperature-2-3.0' u 1:2 w p pt '*'  t 'Ling-model-2',\
'temperature-analy-3.3' u 1:2 w l lw 2 lc rgb 'red' t 'Analytical solution'
pause mouse key