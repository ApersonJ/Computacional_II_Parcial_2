# Examen Parcial II · Física Computacional II (2026-II)

Universidad Distrital Francisco José de Caldas · Programa Académico de Física
Autores: Jesús David Vitola González · Nicolás López Prieto
Profesor: John Hernán Díaz Forero

Este proyecto reúne los tres experimentos del segundo parcial. Todo está escrito en C++ con programación orientada a objetos y se compila y ejecuta desde la raíz con un único `Makefile`.

| Punto | Carpeta | Experimento | Método numérico |
|---|---|---|---|
| 1 · A | `punto1_A_fenomenologia/` | **A2**: reflexión de una partícula contra una pared infinita | Diferencias finitas + Crank–Nicolson (Eigen) |
| 1 · B | `punto1_B_fenomenologia/` | **B1**: oscilador armónico cuántico 1D | Diferencias finitas + Crank–Nicolson (Eigen) |
| 2 | `punto2_motor_molecular/` | **Flashing ratchet**: trinquete que se enciende y se apaga (potencial diente de sierra) | Integrador de Langevin GJF (Grønbech-Jensen–Farago) |

---

## 1. Estructura del proyecto

```
Examen2_FisicaComputacional2_2026II/
├── README.md
├── Makefile                      <- compila y ejecuta los tres puntos
├── punto1_A_fenomenologia/       (A2)
│   ├── include/                  cabeceras .h
│   ├── src/                      .cpp y main.cpp
│   ├── results/                  datos y gráficas (se generan al ejecutar)
│   ├── scripts/                  scripts de Gnuplot
│   └── docs/                     notas y figuras
├── punto1_B_fenomenologia/       (B1, misma organización)
├── punto2_motor_molecular/       (motor molecular, misma organización)
└── informe_general/
    ├── examen2_informe.tex
    ├── examen2_informe.pdf
    └── referencias.bib
```

Al compilar se crean además `build/` (objetos) y `bin/` (ejecutable) dentro de cada punto. Si falta alguna carpeta, `make setup` crea el árbol completo.

### Dónde va cada archivo `.cpp`

Cada punto necesita sus propios archivos, porque A2 y B1 comparten nombres de clase pero tienen implementaciones distintas. Todos los `.cpp` van en `src/` y todas las cabeceras en `include/`. El archivo principal debe llamarse `main.cpp`.

| Carpeta | Archivos de `src/` |
|---|---|
| `punto1_A_fenomenologia` (A2) | `main`, `malla`, `EstadoCuantico`, `Hamiltoniano`, `CrankNicolson`, `Observables`, `Simulacion`, `SolucionAnalitica`, `Convergencia` |
| `punto1_B_fenomenologia` (B1) | los mismos nueve, más `PotencialArmonico` |
| `punto2_motor_molecular` | `main`, `particle`, `motormodel`, `integrator`, `graphics` |

---

## 2. Dependencias

| Herramienta | Para qué | Versión mínima |
|---|---|---|
| `g++` (o `clang++`) | compilar | con soporte C++17 (se usa `std::filesystem`); GCC 9 o superior |
| `make` | automatizar | cualquiera |
| [Eigen](https://eigen.tuxfamily.org) | álgebra lineal del Punto 1 | 3.x |
| `gnuplot` | gráficas y animaciones | 5.x |
| `pdflatex` y `bibtex` | solo para `make informe` | cualquiera |

Versiones con las que se probó este proyecto: `g++ <completar>`, `Eigen <completar>`, `gnuplot <completar>`.

Instalación en Debian/Ubuntu:

```bash
sudo apt install build-essential libeigen3-dev gnuplot
```

Comprobación rápida desde la raíz:

```bash
make check
```

Si Eigen está en otra ruta (por ejemplo en macOS con Homebrew):

```bash
make EIGEN_CFLAGS=-I/opt/homebrew/include/eigen3 all
```

---

## 3. Compilación y ejecución

Todos los comandos se lanzan desde la raíz del proyecto.

| Comando | Qué hace |
|---|---|
| `make help` | muestra la lista de comandos |
| `make check` | verifica compilador, Eigen y gnuplot |
| `make setup` | crea el árbol de carpetas si falta |
| `make all` | compila los tres puntos |
| `make punto1A` · `punto1B` · `punto2` | compila un solo punto |
| `make run1A` · `run1B` | ejecuta el Punto 1 con menú interactivo |
| `make sim1A` · `sim1B` | simulación directa (opción 0) |
| `make conv1A` · `conv1B` | estudio de convergencia (opción 1) |
| `make run2` | ejecuta el Punto 2 y pide los parámetros por teclado |
| `make demo2-det` | Punto 2 determinista con parámetros de ejemplo |
| `make demo2-est` | Punto 2 estocástico con parámetros de ejemplo |
| `make informe` | compila el informe LaTeX a PDF |
| `make clean` | borra `build/` y `bin/` (no toca resultados) |
| `make clean-results` | borra los `.dat`, `.png` y `.gif` generados |

Los comandos de ejecución compilan automáticamente si hace falta. Los de compilación usan `-std=c++17 -O2 -Wall -Wextra`; se pueden cambiar con `make OPT=-O3 all` o `make EXTRA=-g all`.

Cada punto se ejecuta desde su propia carpeta, porque el código escribe en rutas relativas: el Punto 1 en `results/` y `scripts/`, y el Punto 2 en `../results/` y `../scripts/` (se lanza desde `bin/`). El `Makefile` ya se encarga de esto y de crear las carpetas de salida.

---

## 4. Punto 1: Schrödinger 1D (A2 y B1)

Ecuación resuelta, en unidades ħ = m = 1:

```
i ∂ψ/∂t = −(1/2) ∂²ψ/∂x² + V(x) ψ
```

- **Espacio:** diferencias finitas de tres puntos. El Hamiltoniano es tridiagonal, simétrico y disperso (`Eigen::SparseMatrix<double>`).
- **Tiempo:** Crank–Nicolson, `(I + iΔt H/2ħ) ψⁿ⁺¹ = (I − iΔt H/2ħ) ψⁿ`. La matriz de la izquierda se factoriza una sola vez con `Eigen::SparseLU`.
- **Integrales:** sumas con peso Δx.

### Parámetros base

| Parámetro | A2 | B1 |
|---|---|---|
| Potencial | V = 0, ψ(0,t) = 0 | V = ½ω²x² |
| Dominio | [−25, 0] | [−10, 10] |
| Nodos interiores N | 499 | 399 |
| Δx | 0.05 | 0.05 |
| Δt | 0.002 | 0.002 |
| Estado inicial | gaussiana, x₀ = −10, σ = 1, k₀ = 5 | gaussiana, x₀ = 3, σ = √0.5, k₀ = 0 |
| Tiempo final | 4 | 4π |
| Referencia exacta | método de imágenes | estado coherente |

Los parámetros están como constantes al inicio de `src/main.cpp` de cada punto. Para cambiarlos se edita ese archivo y se vuelve a compilar con `make`.

### Menú del programa

Al ejecutarlo pide un número:

- `0`: simulación completa (series, densidades, Ehrenfest y gráficas).
- `1`: estudio de convergencia en Δx y Δt, y comparación con la solución exacta.

### Archivos de salida

Todos en `results/` del punto correspondiente.

| Archivo | Columnas |
|---|---|
| `series.dat` | `t`, norma, ⟨x⟩, ⟨pₓ⟩, ⟨H⟩, probabilidad de una región (A2: P(−5<x<0) · B1: P(x>0)) |
| `densidad.dat` | `x`, `t`, ρ (un bloque por instante, separados por línea en blanco) |
| `ehrenfest.dat` | `t`, ⟨pₓ⟩, m d⟨x⟩/dt, diferencia |
| `convergencia/espacial.dat` | Δx, error L² de ρ, \|⟨H⟩ − ⟨H⟩ exacta\|, deriva de ⟨H⟩, \|N − 1\| |
| `convergencia/temporal.dat` | Δt, error L² de ρ frente a la corrida de referencia, deriva de ⟨H⟩, \|N − 1\| |
| `convergencia/comparacion_densidad.dat` | (solo A2) `x`, `t`, ρ numérica, ρ analítica |
| `convergencia/fit.log` | ajuste log-log (pendientes) de Gnuplot |
| `graficas/*.png` | gráficas generadas por Gnuplot |

Los scripts de Gnuplot se escriben en `scripts/` (por ejemplo `graficar.plt` y `graficar_convergencia.plt`), así que las gráficas se pueden regenerar a mano con `gnuplot scripts/graficar.plt` desde la carpeta del punto.

### Clases

| Clase | Responsabilidad |
|---|---|
| `Malla` | nodos interiores y Δx |
| `EstadoCuantico` | amplitudes complejas ψⱼ (`VectorXcd`), gaussiana, norma |
| `Hamiltoniano` | ensambla H dispersa y verifica simetría |
| `CrankNicolson` | matrices A y B, factorización y avance temporal |
| `Observables` | densidad, ⟨x⟩, ⟨pₓ⟩, ⟨H⟩ y probabilidad por región |
| `Simulacion` | ciclo temporal, archivos de salida y script de Gnuplot |
| `SolucionAnalitica` | referencia exacta (imágenes en A2, coherente en B1) |
| `Convergencia` | refinamiento en Δx y Δt |
| `PotencialArmonico` | (solo B1) V(x) = ½ mω²x² evaluado en los nodos |

---

## 5. Punto 2: flashing ratchet

Partícula en un potencial de diente de sierra periódico que se enciende y se apaga en el tiempo:

```
U(x,t) = m(t) · U_sierra(x)        m(t) = 1 durante τ_on, 0 durante τ_off
```

El potencial tiene periodo L, profundidad U₀ y asimetría α. Se integra la ecuación de Langevin con el esquema GJF, con unidades kB = 1 y paso Δt = 0.001. Se guardan datos cada 10 pasos.

### Parámetros

El programa los pide por teclado, en este orden:

| # | Parámetro |
|---|---|
| 1 | Modelo: `1` estocástico (fricción y ruido) · `2` determinista |
| 2 | Duración de la simulación |
| 3 | Duración del potencial encendido (τ_on) |
| 4 | Duración del potencial apagado (τ_off) |
| 5 | Posición inicial x₀ |
| 6 | Velocidad inicial v₀ |
| 7 | Masa m |
| 8 | Periodo del potencial L |
| 9 | Profundidad del pozo U₀ |
| 10 | Parámetro de asimetría α (entre 0 y 0.5) |
| 11 | *(solo estocástico)* coeficiente de fricción γ |
| 12 | *(solo estocástico)* energía térmica kB·T |

En el modelo estocástico se promedian 500 realizaciones del ruido (semillas 1000 + k). En el determinista hay una sola trayectoria. Los valores de `demo2-det` y `demo2-est` son solo un ejemplo; se pueden cambiar en las variables `DEMO_DET` y `DEMO_EST` del `Makefile`.

### Archivos de salida

En `punto2_motor_molecular/results/`:

| Archivo | Columnas |
|---|---|
| `trayectoria_individual.dat` | `t`, `x`, `v`, `p`, `K`, `U`, `E` (primera realización) |
| `promedios_ensamble.dat` | `t`, ⟨x⟩, ⟨v⟩, ⟨p⟩, ⟨K⟩, ⟨U⟩, ⟨E⟩ (promedio del ensamble) |
| `animations/*.gif` | `position`, `phasespace`, `velocity`, `momentum`, `totalenergy`, `meanposition`, `meankinetic`, `meanpotential` |

Al terminar, el programa imprime la velocidad media estacionaria ⟨v⟩ = (x_final − x₀)/T promediada sobre el ensamble. Los scripts de Gnuplot de las animaciones (`position.gp`, etc.) quedan en `scripts/`.

### Clases

| Clase | Responsabilidad |
|---|---|
| `Particle` | posición, velocidad, masa, momentum |
| `MotorModel` | potencial de diente de sierra, fuerza, energías, fricción y temperatura |
| `Integrator` | paso GJF y generación del ruido |
| `Graphics` | scripts y animaciones de Gnuplot |

---

## 6. Cómo reproducir todos los resultados

```bash
# 1. Comprobar dependencias y compilar todo
make check
make all

# 2. Punto 1-A (A2)
make sim1A          # simulación
make conv1A         # convergencia y comparación con la solución exacta

# 3. Punto 1-B (B1)
make sim1B
make conv1B

# 4. Punto 2
make demo2-det      # determinista, parámetros de ejemplo
make demo2-est      # estocástico, parámetros de ejemplo
# o bien, con los parámetros del informe:
make run2

# 5. Informe
make informe
```

Los resultados quedan en la carpeta `results/` de cada punto. Para empezar desde cero: `make clean clean-results`.

---

## 7. Notas

- La conversión de gráficas depende de `gnuplot`. Si no está instalado, los programas terminan igual, avisan del error y dejan los `.dat` y los scripts listos para graficar después.
- El tiempo de ejecución depende de Δx y Δt (Punto 1) y de las 500 realizaciones (Punto 2 estocástico).
- La solución de referencia de A2 y B1 solo valida el cálculo: no sustituye a Crank–Nicolson.
