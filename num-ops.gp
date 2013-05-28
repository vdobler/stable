#!/usr/bin/gnuplot -persist
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
set timefmt z "%d/%m/%y,%H:%M"
set zdata 
set timefmt y "%d/%m/%y,%H:%M"
set ydata 
set timefmt x "%d/%m/%y,%H:%M"
set xdata 
set timefmt cb "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set y2data 
set timefmt x2 "%d/%m/%y,%H:%M"
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy x,y
set format x "%.0e"
set format y "%.0e"
set format x2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
set grid nopolar
set grid xtics nomxtics ytics nomytics noztics nomztics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000
set raxis
set key title ""
set key inside left top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1.5 width 0 height 1
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
unset logscale
set logscale x 10
set logscale y 10
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics 10.000000
set mytics default
set mztics default
set mx2tics default
set my2tics 0.000000
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
# set y2tics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
# set y2tics autofreq  norangelimit
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "{/*1.3 Stable Sorting vs. Standard Sort }\nMergeSort (SymMerge/BlockSwap) vs. QuickSort" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "{/*1.25 n}" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ 100.000 : 1.00000e+09 ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "{/*1.25 # Swaps}" 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label ""
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ 1e2 : 1e12 ] noreverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_US.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit noerrorvariables
GNUTERM = "wxt"
a = 1.46773740921261
FIT_CONVERGED = 1
FIT_NDF = 11
FIT_STDFIT = 71.1228959183182
FIT_WSSR = 55643.1295618872
b = 7.48996227352459
c = 1.16239563551843
d = 1.001

# m = Merge Sort (Stable)
fit [1e4:] mmb*(x**mmc) "num-ops" using 1:2:($2/1000) via mmb,mmc
fit [1e4:] mmal*x*log(x) "num-ops" using 1:2:($2/1000) via mmal
fit [1e4:] mma*x*log(x) "num-ops" using 1:3:($3/1000) via mma
fit [1e4:] mub*(x**muc) "num-ops" using 1:4:($4/1000) via mub,muc
fit [1e4:] mual*x*log(x) "num-ops" using 1:4:($4/1000) via mual
fit [1e4:] mua*x*log(x) "num-ops" using 1:5:($5/1000) via mua

# q = Quick Sort (Sort)
fit [1e4:] qmal*x*log(x) "num-ops" using 1:6:($6/1000) via qmal
fit [1e4:] qma*x*log(x) "num-ops" using 1:7:($7/1000) via qma
fit [1e4:] qual*x*log(x) "num-ops" using 1:8:($8/1000) via qual
fit [1e4:] qua*x*log(x) "num-ops" using 1:9:($9/1000) via qua

# k = Katajainen
fit ka*x*log(x) "num-ops" using 1:10:($10/1000) via ka
fit kl*x*log(x) "num-ops" using 1:11:($10/1000) via kl

tmms = sprintf("Stable Multiple: %.1f n^{/*1.1 %.2f}", mmb,mmc)
tmml = sprintf("Stable Multiple: %.2f n log n", mma)
tmus = sprintf("Stable Unique: %.1f n^{/*1.1 %.2f}", mub,muc)
tmul = sprintf("Stable Unique: %.2f n log n", mua)

tqms = sprintf("Sort Multiple: %.2f n log n", qmal)
tqml = sprintf("Sort Multiple: %.2f n log n", qma)
tqus = sprintf("Sort Unique: %.2f n log n", qual)
tqul = sprintf("Sort Unique: %.2f n log n", qua)

set term pngcairo enhanced color solid font "Helvetica" fontscale 0.95 size 1000,800

set label "Distributions:" at screen 0.75,0.22
set label "Unique = rand.Intn(n * 10)" at screen 0.78,0.2
set label "Multiple = rand.Intn(n / 100)" at screen 0.777,0.18

set out "num-swaps.png"
plot "num-ops" u 1:2 notit with points ls 1,  mmb*(x**mmc) tit tmms with lines ls 1, \
     "" u 1:4 notit with points ls 2, mub*(x**muc) tit tmus with lines ls 2,         \
     "" u 1:6 notit with points ls 3, qmal*x*log(x) tit tqms with lines ls 3,         \
     "" u 1:8 notit with points ls 7, qual*x*log(x) tit tqus with lines ls 7, \
     mmal*x*log(x) tit "n log n - fit to Stable", \
     "" u 1:10 tit "Katajainen 4*-mergesort" with points ls 4, \
     ka*x*log(x) tit "Katajainen n log n fit" with lines ls 4


set ylabel "{/*1.25 # Less}" 
set yrange [ 1e3 : 1e11 ]
set mytics 0
set out "num-less.png"
plot "num-ops" u 1:3 notit with points ls 1, mma*x*log(x) tit tmml with lines ls 1, \
     "" u 1:5 notit with points ls 2, mua*x*log(x) tit tmul with lines ls 2, \
     "" u 1:7 notit with points ls 3, qma*x*log(x) tit tqml with lines ls 3, \
     "" u 1:9 notit with points ls 7, qua*x*log(x) tit tqul with lines ls 7, \
     "" u 1:11 tit "Katajainen 4*-mergesort" with points ls 4

unset logscale y
set yrange [0:30]
set format y "%.0f"
set mytics 5
set ylabel "{/*1.25 (# Stable) / (# Sort)}"
set out "num-ratio.png"
plot "num-ops" u 1:($2/$6) tit "Swaps Multiple",  \
     "" u 1:($3/$7) tit "Less Multiple",          \
     "" u 1:($4/$8) tit "Swaps Unique",           \
     "" u 1:($5/$9) tit "Less Unique"


