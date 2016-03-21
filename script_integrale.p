set xrange [2000:1]
set logscale x
set logscale y

plot 'energy_rate.dat' u 1:2 w l , \
'energy_rate_2.dat' u 1:2 w l
