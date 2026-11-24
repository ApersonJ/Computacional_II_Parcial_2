datafile = '../results/promedios_ensamble.dat'
set terminal gif animate delay 4 loop 0 size 900,500
set output '../results/animations/meanposition.gif'
# tiempo (columna 1) para el titulo de cada cuadro
stats datafile using 1 nooutput
N    = STATS_records
tmin = STATS_min
tmax = STATS_max
# rangos fijos para que los ejes no cambien entre cuadros
stats datafile using 1:2 nooutput
xmin = STATS_min_x
xmax = STATS_max_x
ymin = STATS_min_y
ymax = STATS_max_y
padx = 0.05*(xmax - xmin)
pady = 0.05*(ymax - ymin)
if (padx == 0) padx = 1
if (pady == 0) pady = 1

set xrange [xmin-padx:xmax+padx]
set yrange [ymin-pady:ymax+pady]
set xlabel 't'
set ylabel '<x>(t)'
set grid
set key top left

# maximo 150 cuadros
nframes = (N-1 > 150) ? 150 : N-1

do for [k=0:nframes] {
    idx = (k*(N-1))/nframes
    ti = tmin + (tmax-tmin)*idx/(N-1)
    set title sprintf('Posicion media del ensamble  (t = %.2f)', ti)
    plot datafile every ::0::idx using 1:2 with lines lw 2 lc rgb '#1f77b4' title '<x>(t)', datafile every ::idx::idx using 1:2 with points pt 7 ps 1.8 lc rgb '#d62728' title 'actual'
}

unset output
