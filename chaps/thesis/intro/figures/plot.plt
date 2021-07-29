#! /usr/bin/env gnuplot

set term pngcairo size 1000,400 font ",14"
set output "distrib.png"

Theads = 'heads size 0.2,90 front ls 201'

M(x,r,u0,T) = r/(sqrt(2*pi*T))*exp(-(x-u0)**2/(2*T))

a  = 0.2
Tc = 0.2
v0 = 5
fc(x) = M(x,1-a,0,Tc)
fh(x) = M(x,a/2,v0,1) + M(x,a/2,-v0,1)

x0 = 0.6
set arrow from -x0+0.05,fc(-x0) to x0-0.05,fc(x0) @Theads lw 2 lc rgb "#242424"
set label 1 at 0-0.07,fc(x0)+0.05 "T_c" front

set label 2 at -v0-0.15,fh(-v0)+0.1 "-v_0" front
set label 3 at  v0-0.15,fh(v0) +0.1 " v_0" front

set samples 1000
set xlabel "v"

plot fc(x)+fh(x) w l lw 2 lc rgb "#535c68" notitle

