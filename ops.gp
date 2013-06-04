set title "Sort vs Stable"

set logscale x
set xlabel "N"
set xrange [3e0:3e9]
set format x "%.0e"

set logscale y
set ylabel "# Swap"
set yrange [2e0:1e12]
set format y "%.0e"

set key left

set term dumb 135,43
set mxtics 0; set mytics 0

set ytics 10 nomirror
set mxtics 10
set mytics 10

set grid


nlogn(a,n) = a * n * log(n)
nlog2n(a,n) = a * n * log(n) * log(n)


fit nlog2n(StableSwap,x) "ops" u 1:3:(sqrt($3)) via StableSwap
fit nlogn(SortSwap,x) "ops" u 1:7:(sqrt($7)) via SortSwap

fit nlogn(StableLess,x) "ops" u 1:5:(sqrt($5)) via StableLess
fit nlogn(SortLess,x) "ops" u 1:9:(sqrt($5)) via SortLess


plot nlog2n(StableSwap,x) notit ls 1,            \
     nlogn(SortSwap,x) notit ls 2,               \
     "ops" u 1:3 tit "Stable: Swap" ls 1 ps 2,   \
     "ops" u 1:7 tit "Sort: Swap" ls 2 ps 1.5,   \
     nlogn(StableLess,x)  notit ls 3,            \
     nlogn(SortLess,x) notit ls 4,               \
     "ops" u 1:5 tit "Stable: Less" ls 3 ps 1.5, \
     "ops" u 1:9 tit "Sort: Less" ls 4 ps 1.5

set term pngcairo enhanced color solid font "Helvetica" fontscale 0.95 size 1000,800
set out "ops.png"
replot

