# set terminal png
# set view map
# set rmargin at screen 0.8
# do for [i=1:10]{
#   set output sprintf('Density_n%d.png', i)
#   sp sprintf('Density_n%d.dat', i) u 1:2:($3 + $4 + $5 + $6 + $7 + $8) palette z pt 7 ps 3
# }


# reset
# set terminal png
# set view map
# set rmargin at screen 0.8
# do for [i=1:10]{
#   set output sprintf('Psi_1_n%d.png', i)
#   sp sprintf('Psi_1_n%d.dat', i) u 1:2:($3 + $4 + $5 + $6 + $7 + $8) palette z pt 7 ps 3
# }

reset
set terminal png
set output 'E_1(B).png'
set xlabel 'B (T)'
set ylabel 'E_1 (meV)'
set title 'ℏ{/Symbol w} = 37.378 meV'
plot for [i = 0:5] sprintf('Energies1_B%d.dat', i*10) u (10*i):2 pt 7 ps 1 lc rgb 'black' notitle

reset
set terminal png
set output 'E_2(B).png'
set xlabel 'B (T)'
set ylabel 'E_2 (meV)'
set title 'ℏ{/Symbol w} = 37.378 meV'
plot for [i = 0:5] sprintf('Energies2_B%d.dat', i*10) u (10*i):2 pt 7 ps 1 lc rgb 'black' notitle


reset
set terminal png
set output 'E_1_om18(B).png'
set xlabel 'B (T)'
set ylabel 'E_1 (meV)'
set title 'ℏ{/Symbol w} = 18.689 meV'
plot for [i = 0:5] sprintf('Energies1_B%d_om18.dat', i*10) u (10*i):2 pt 7 ps 1 lc rgb 'black' notitle

reset
set terminal png
set output 'E_2_om18(B).png'
set xlabel 'B (T)'
set ylabel 'E_2 (meV)'
set title 'ℏ{/Symbol w} = 18.689 meV'
plot for [i = 0:5] sprintf('Energies2_B%d_om18.dat', i*10) u (10*i):2 pt 7 ps 1 lc rgb 'black' notitle