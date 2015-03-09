set logscale
set xrange [20:1000]
plot 'energy_rate_slatyer_avec_halo_ee1GeV.dat' u 1:3, \
'energy_rate_vivian_avec_halo_ee1GeV.dat' u 1:3 , \
'energy_rate_slatyer_avec_halo_ee1TeV.dat' u 1:3, \
'energy_rate_vivian_avec_halo_ee1TeV.dat' u 1:3
