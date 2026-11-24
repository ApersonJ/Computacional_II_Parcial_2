# Generado por Simulacion::EscribirScriptGnuplot()
dt_foto = 0.25
n_fotos = 17

set terminal pngcairo size 900,600 enhanced
set grid
set key top right

# --- Series temporales ---
set xlabel 't'

# Norma: casi constante en 1, rango muy estrecho y con decimales visibles.
set output 'results/graficas/norma.png'
set ylabel 'N(t)'
set yrange [1-1e-6:1+1e-6]
set format y '%.8f'
plot 'results/series.dat' using 1:2 with lines lw 2 title 'norma'
set format y '% h'
set autoscale y

set output 'results/graficas/x_medio.png'
set ylabel '<x>(t)'
plot 'results/series.dat' using 1:3 with lines lw 2 title '<x>'

set output 'results/graficas/px_medio.png'
set ylabel '<p_x>(t)'
plot 'results/series.dat' using 1:4 with lines lw 2 title '<p_x>'

# Energía: casi constante, se deja el rango automático.
set output 'results/graficas/energia.png'
set ylabel '<H>(t)'
set format y '%.6f'
plot 'results/series.dat' using 1:5 with lines lw 2 title '<H>'
set format y '% h'

set output 'results/graficas/prob_cerca_pared.png'
set ylabel 'P(-5 < x < 0)'
plot 'results/series.dat' using 1:6 with lines lw 2 title 'prob. cerca de la pared'

# --- Perfiles de densidad en varios instantes ---
set output 'results/graficas/densidad_perfiles.png'
set xlabel 'x'
set ylabel '{/Symbol r}(x,t) = |{/Symbol y}|^2'
plot for [i=0:n_fotos-1:4] 'results/densidad.dat' every :::i::i using 1:3 with lines lw 2 title sprintf('t = %.2f', i*dt_foto)

# --- Ehrenfest: <px> frente a m d<x>/dt ---
set output 'results/graficas/ehrenfest.png'
set xlabel 't'
set ylabel '<p_x>  y  m d<x>/dt'
plot 'results/ehrenfest.dat' using 1:2 with lines lw 3 title '<p_x> (derivada espacial)', \
     'results/ehrenfest.dat' using 1:3 with points pt 6 ps 0.8 title 'm d<x>/dt (derivada temporal)'

set output 'results/graficas/ehrenfest_diferencia.png'
set ylabel '<p_x> - m d<x>/dt'
plot 'results/ehrenfest.dat' using 1:4 with lines lw 2 title 'diferencia'

# --- Mapa de calor de rho(x,t) ---
set output 'results/graficas/mapa_densidad.png'
set xlabel 'x'
set ylabel 't'
set pm3d map
set palette defined (0 'white', 1 'blue', 2 'red')
splot 'results/densidad.dat' using 1:2:3 notitle
