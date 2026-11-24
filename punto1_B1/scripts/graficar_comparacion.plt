tiempos = "1 2 3 "

set terminal pngcairo size 900,600 enhanced
set grid
set key top left
set xlabel 'x'
set ylabel '{/Symbol r}(x,t) = |{/Symbol y}|^2'
set xrange [-15:0]

do for [i=0:words(tiempos)-1] {
  t = word(tiempos, i+1)
  set output sprintf('results/graficas/comparacion_densidad_t%s.png', t)
  set title sprintf('Densidad numerica y analitica en t = %s', t)
  plot 'results/convergencia/comparacion_densidad.dat' every :::i::i using 1:3 with lines linewidth 2 title 'Crank-Nicolson', \
       '' every 4:::i::i using 1:4 with points pointtype 6 pointsize 0.7 title 'metodo de imagenes'
}
