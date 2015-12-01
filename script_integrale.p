set xrange [2000:1]
set yrange [0.5:1]
set logscale x
set logscale y

plot 'integrale_10000.dat' u 1:2 w l , \
'integrale_1000.dat' u 1:2 w l , \
'integrale_100.dat' u 1:2 w l , \
'integrale_10.dat' u 1:2 w l, \
'integrale_hyrec.dat' u 1:2 
