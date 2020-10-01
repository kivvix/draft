#! /usr/bin/env gnuplot

set term png size 2400,800 font "arial,22."
set output "fend.png"

set xrange [-5:5]
set yrange [-5:5]
set zrange [0:1]
set cbrange [0:1]

set xlabel "x"
set ylabel "y"

set multiplot layout 1,3 \
    margins 0.05,0.93,0.15,0.85 \
    spacing 0.03,0.03

set title "f"
unset key
unset colorbox
set ytics
plot "fend.dat" u 1:2:3 w image

set title "g"
unset key
unset colorbox
unset ylabel; unset ytics
plot "gend.dat" u 1:2:3 w image

set title "g filtré (donc f...)"
unset key
set colorbox
unset ylabel; unset ytics
plot "ffend.dat" u 1:2:3 w image

set colorbox user origin 0.15,0.1 size 0.05,0.8
unset multiplot

set term png size 800,800 font "arial,22."
set output "finit.png"
set ytics
set ylabel
unset colorbox
set title "f(t=0)"
plot "finit.dat" u 1:2:3 w image


