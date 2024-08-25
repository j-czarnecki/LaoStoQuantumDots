set terminal png
set output 'density.png'

set view map
set rmargin at screen 0.8
sp 'Density_n1.dat' u 1:2:($3 + $4 + $5 + $6 + $7 + $8) palette z pt 7 ps 2


reset
set terminal png
set output 'Psi_1.png'

set view map
set rmargin at screen 0.8
sp 'Psi_1_n4.dat' u 1:2:($3+$4) palette z pt 7 ps 2