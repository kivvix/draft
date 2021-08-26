#! /usr/bin/env gnuplot

#set term pngcairo size 1000,400 font ",14"
#set output "distrib.png"
set term epslatex color font ",5" header \
   "\\newcommand{\\ft}[0]{\\footnotesize}"
set size 1.125,0.45
set output "distrib.tex"

Theads = 'heads size 0.2,90 front ls 201'

M(x,r,u0,T) = r/(sqrt(2*pi*T))*exp(-(x-u0)**2/(2*T))

a  = 0.2
Tc = 0.2
vo = 5
fc(x) = M(x,1-a,0,Tc)
fh(x) = M(x,a/2,vo,1) + M(x,a/2,-vo,1)

xo = 0.6
set arrow from -xo+0.05,fc(-xo) to xo-0.05,fc(xo) @Theads lw 2 lc rgb "#242424"
set label 1 at 0-0.5,fc(xo)+0.07 "\\ft $\\sqrt{T_c}$" front font ",3"

set label 2 at -vo-0.6,fh(-vo)+0.1 "$-v_0$" front
set label 3 at  vo-0.3,fh(vo) +0.1 "$v_0$"  front

set samples 1000
set xlabel "$v$" offset 0,graph -0.1

set xtics offset 0,graph -0.1
set ytics 0.25

plot fc(x)+fh(x) w l lw 3 lc rgb "#535c68" notitle

set output # Closes the temporary output files.
!gsed -i 's/distrib/\\localPath\/figures\/distrib/g' distrib.tex
