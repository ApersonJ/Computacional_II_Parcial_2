# Generado por Convergencia::Graficar()
T = 2

set terminal pngcairo size 900,600 enhanced
set grid
set logscale xy
set format xy '10^{%L}'
set fit quiet
set fit logfile 'results/convergencia/fit.log'

# Ajuste lineal en log-log: log10(error) = p * log10(h) + b
f_esp(x) = p_esp*x + b_esp
f_ene(x) = p_ene*x + b_ene
f_tem(x) = p_tem*x + b_tem
fit f_esp(x) 'results/convergencia/espacial.dat' using (log10($1)):(log10($2)) via p_esp,b_esp
fit f_ene(x) 'results/convergencia/espacial.dat' using (log10($1)):(log10($3)) via p_ene,b_ene
fit f_tem(x) 'results/convergencia/temporal.dat' using (log10($1)):(log10($2)) via p_tem,b_tem

set output 'results/graficas/convergencia_espacial.png'
set xlabel '{/Symbol D}x'
set ylabel 'error L2 de {/Symbol r}(x,T) frente a la solucion analitica'
set title sprintf('Refinamiento espacial en t = %.1f', T)
plot 'results/convergencia/espacial.dat' using 1:2 with linespoints pointtype 7 linewidth 2 \
     title sprintf('numerico (pendiente = %.2f)', p_esp)

set output 'results/graficas/convergencia_energia.png'
set ylabel '|<H> - <H>_{continuo}|'
set title 'Error de la energia media frente a {/Symbol D}x'
plot 'results/convergencia/espacial.dat' using 1:3 with linespoints pointtype 7 linewidth 2 \
     title sprintf('numerico (pendiente = %.2f)', p_ene)

set output 'results/graficas/convergencia_temporal.png'
set xlabel '{/Symbol D}t'
set ylabel 'error L2 de {/Symbol r}(x,T) frente a la corrida de referencia'
set title sprintf('Refinamiento temporal en t = %.1f', T)
plot 'results/convergencia/temporal.dat' using 1:2 with linespoints pointtype 7 linewidth 2 \
     title sprintf('numerico (pendiente = %.2f)', p_tem)
