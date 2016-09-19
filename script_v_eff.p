set xrange [1:2000]
set logscale x
set logscale y
f1(x)=b*x**a
f2(x)=c*x**d
fit [x = 10:180]  f1(x) "v_eff.dat" using 1:2 via a,b
fit [x = 180:1500]  f2(x) "v_eff.dat" using 1:2 via c,d
plot 'v_eff.dat' u 1:2 w l, \
f1(x), \
f2(x)
