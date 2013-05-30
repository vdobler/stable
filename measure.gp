# set term dumb 120,40
set term png size 1200,840
set xrange [0:40]
set xtics 4

set out "measure-2n.png"
set multiplot layout 2,2
n = 100
while (n < 1e6) {
  set title sprintf("2*N : n = %d", n);
  k = n/100
  plot "2N" u 2:($1==n ? $3/k : 1/0) tit "Rec", "2N" u 2:($1==n ? $4/k : 1/0) tit "BU"
  n = n *10;
}
unset multiplot

set out "measure-n2.png"
set multiplot layout 2,2
n = 100
while (n < 1e6) {
  set title sprintf("N/2 : n = %d", n);
  k = n/100
  plot        "N2" u 2:($1==n ? $3/k : 1/0) tit "Rec", "N2" u 2:($1==n ? $4/k : 1/0) tit "BU"
  n = n *10;
}
unset multiplot



