# Guardar imagen PNG
set terminal pngcairo size 1200,900 enhanced font 'Verdana,10'
set output 'resultado_gnuplot.png'

set title 'Distribución de Voltaje'
set xlabel 'X (m)'
set ylabel 'Y (m)'
set zlabel 'Voltaje (V)'

set grid
set dgrid3d 50,50 splines
set hidden3d
set view 60, 60
set contour base
set cntrparam levels incremental 0,10,100
set style data lines
set datafile commentschars '#'

splot 'laplace.dat' using 1:2:3 with lines lc rgb 'black' title 'V(x,y)', \
          '' using 1:2:3 with pm3d notitle 

# Mostrar en pantalla usando wxt (si hay entorno gráfico)
	      set terminal wxt enhanced
	      unset output
	      replot
