set logscale
set xrange [1:4000]
plot 'energy_rate_1Msun.dat' u 5:3, \
'energy_rate_10Msun.dat' u 5:3, \
'energy_rate_100Msun.dat' u 5:3, \
'energy_rate_1000Msun.dat' u 5:3
