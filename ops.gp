set title "{/*1.3 Sort vs. Stable}"

set logscale x
set xlabel "N"
set xrange [1e1:1e8]
set format x "%.0e"

set logscale y
set ylabel "Calls to Swap and Less"
set yrange [1e1:1e11]
set format y "%.0e"

unset logscale y2
set y2label "Swaps_{Stable} / Swaps_{Sort}"
set y2range [0:25]
set y2tics 5
set my2tics 5
set format y2 "%g"


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
     "ops" u 1:9 tit "Sort: Less" ls 4 ps 1.5,   \
     "ops" u 1:($3/$7) ax x1y2 tit "Swap ratio" ls 5

set term pngcairo enhanced color solid font "Helvetica" fontscale 0.95 size 1000,800
set out "ops.png"
replot

