set terminal png

set view map
set rmargin at screen 0.8
do for [i=1:10]{
  set output sprintf('Density_n%d.png', i)
  sp sprintf('Density_n%d.dat', i) u 1:2:($3 + $4 + $5 + $6 + $7 + $8) palette z pt 7 ps 3
}


reset
set terminal png
set view map
set rmargin at screen 0.8
do for [i=1:10]{
  set output sprintf('Psi_1_n%d.png', i)
  sp sprintf('Psi_1_n%d.dat', i) u 1:2:($3 + $4 + $5 + $6 + $7 + $8) palette z pt 7 ps 3
}