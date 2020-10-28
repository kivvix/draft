#! /usr/bin/env gnuplot

func="cos sin sinHF hat gat"
scheme="cd2 weno"

plot for [j=1:words(func)] for [i=1:words(scheme)] "end_".word(scheme,i)."_".word(func,j).".dat" u 1:3 w p title word(func,j)." ".word(scheme,i)
pause -1
plot for [j=1:words(func)] for [i=1:words(scheme)] "end_".word(scheme,i)."_".word(func,j).".dat" u 1:2 w lp title word(func,j)." ".word(scheme,i)
pause -1
