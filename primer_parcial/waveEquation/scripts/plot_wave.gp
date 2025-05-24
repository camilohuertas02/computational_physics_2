# script.gnuplot

# Parámetro de tiempo
t = 0.123   # <-- Cambia este valor según el tiempo que desees graficar

set title sprintf("Onda analítica vs numérica: t = %.3f s", t)
set xlabel 'x'
set ylabel 'y(x,t)'
set xrange [0:4]
set yrange [-2.5:2.5]
set grid

plot 'dataA.dat' using 1:2 with lines lw 2 lc rgb 'blue' dt 2 title 'Analítica', \
     'dataN.dat' using 1:2 with lines lw 2 lc rgb 'red' title 'Numérica'
