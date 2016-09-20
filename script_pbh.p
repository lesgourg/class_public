set logscale
set xrange [1:4000]
plot 'energy_rate_PBH1Msun_v2.dat' u 7:2, \
'energy_rate_PBH10Msun_v2.dat' u 7:2, \
'energy_rate_PBH100Msun_v2.dat' u 7:2, \
'energy_rate_PBH1e3Msun_v2.dat' u 7:2, \
'energy_rate_PBH1e4Msun_v2.dat' u 7:2, \
'energy_rate_PBH1e5Msun_v2.dat' u 7:2, \
'energy_rate_PBH1e6Msun_v2.dat' u 7:2, \
'energy_rate_PBH1e7Msun_v2.dat' u 7:2, \
'energy_rate_PBH1e8Msun_v2.dat' u 7:2
