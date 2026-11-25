# =====================================================================
# Makefile general - Examen Parcial II, Fisica Computacional II (2026-II)
# Vive en la raiz del proyecto, fuera de las carpetas de los puntos.
#   make help   -> lista de comandos
# =====================================================================
.DEFAULT_GOAL := help

CXX      ?= g++
STD      := -std=c++17
OPT      ?= -O2
WARN     ?= -Wall -Wextra
EXTRA    ?=
CXXFLAGS  = $(STD) $(OPT) $(WARN) $(EXTRA)

# Eigen: se busca con pkg-config; si no esta, se usa la ruta tipica de Linux.
# Otra ruta:  make EIGEN_CFLAGS=-I/ruta/a/eigen3
EIGEN_CFLAGS ?= $(shell pkg-config --cflags eigen3 2>/dev/null)
ifeq ($(strip $(EIGEN_CFLAGS)),)
EIGEN_CFLAGS := -I/usr/include/eigen3
endif

PA := punto1_A_fenomenologia
PB := punto1_B_fenomenologia
P2 := punto2_motor_molecular
INF := informe_general

BIN_A := $(PA)/bin/punto1A
BIN_B := $(PB)/bin/punto1B
BIN_2 := $(P2)/bin/motor

# Valores de EJEMPLO para las demos del Punto 2 (ajustalos a los de tu informe).
# Orden de lectura del programa: eleccion, duracion, tau_on, tau_off, x0, v0,
# masa, L, U0, alpha  (+ gamma y energia termica si es estocastico).
DEMO_DET := 2 50 1 1 0.2 0 1 1 5 0.2
DEMO_EST := 1 50 1 1 0.2 0 1 1 5 0.2 1 0.5

# ---------------------------------------------------------------------
# Plantilla de compilacion: $(1)=carpeta  $(2)=ejecutable  $(3)=flags extra
# ---------------------------------------------------------------------
define BUILD_RULES
$(1)_SRCS := $$(wildcard $(1)/src/*.cpp)
$(1)_OBJS := $$(patsubst $(1)/src/%.cpp,$(1)/build/%.o,$$($(1)_SRCS))
$(2): $$($(1)_OBJS) | $(1)/bin
	$$(CXX) $$(CXXFLAGS) $$^ -o $$@
$(1)/build/%.o: $(1)/src/%.cpp | $(1)/build
	$$(CXX) $$(CXXFLAGS) -I$(1)/include -I$(1)/src $(3) -MMD -MP -c $$< -o $$@
$(1)/bin $(1)/build:
	mkdir -p $$@
-include $$($(1)_OBJS:.o=.d)
endef

$(eval $(call BUILD_RULES,$(PA),$(BIN_A),$(EIGEN_CFLAGS)))
$(eval $(call BUILD_RULES,$(PB),$(BIN_B),$(EIGEN_CFLAGS)))
$(eval $(call BUILD_RULES,$(P2),$(BIN_2),))

.PHONY: help check setup all punto1A punto1B punto2 \
        run1A sim1A conv1A run1B sim1B conv1B run2 demo2-det demo2-est \
        informe clean clean-results

help:
	@echo "Uso: make <comando>"
	@echo ""
	@echo "  Preparacion"
	@echo "    check        verifica g++, Eigen y gnuplot"
	@echo "    setup        crea el arbol de carpetas del examen (si faltan)"
	@echo ""
	@echo "  Compilar"
	@echo "    all          compila los tres puntos"
	@echo "    punto1A      compila Punto 1-A (A2: pared infinita)"
	@echo "    punto1B      compila Punto 1-B (B1: oscilador armonico)"
	@echo "    punto2       compila Punto 2 (motor molecular)"
	@echo ""
	@echo "  Ejecutar (compilan si hace falta)"
	@echo "    run1A / run1B   menu interactivo: 0 = simulacion, 1 = convergencia"
	@echo "    sim1A / sim1B   simulacion directa (opcion 0)"
	@echo "    conv1A / conv1B estudio de convergencia (opcion 1)"
	@echo "    run2            Punto 2 interactivo (pide los parametros)"
	@echo "    demo2-det       Punto 2 determinista con parametros de ejemplo"
	@echo "    demo2-est       Punto 2 estocastico con parametros de ejemplo"
	@echo ""
	@echo "  Otros"
	@echo "    informe         compila el informe LaTeX (informe_general/)"
	@echo "    clean           borra build/ y bin/ de los tres puntos"
	@echo "    clean-results   borra .dat, .png y .gif generados en results/"

check:
	@command -v $(CXX) >/dev/null 2>&1 && echo "[ok]    $(CXX)" || echo "[FALTA] $(CXX)"
	@command -v gnuplot >/dev/null 2>&1 && echo "[ok]    gnuplot" || echo "[FALTA] gnuplot (las graficas no se generaran)"
	@echo '#include <Eigen/Dense>' | $(CXX) $(EIGEN_CFLAGS) -x c++ -fsyntax-only - >/dev/null 2>&1 && echo "[ok]    Eigen ($(EIGEN_CFLAGS))" || echo "[FALTA] Eigen: instalalo o usa make EIGEN_CFLAGS=-I/ruta"
	@command -v pdflatex >/dev/null 2>&1 && echo "[ok]    pdflatex" || echo "[aviso] pdflatex no encontrado (solo para make informe)"

setup:
	for p in $(PA) $(PB) $(P2); do mkdir -p $$p/include $$p/src $$p/results $$p/scripts $$p/docs; done
	mkdir -p $(INF)
	@echo "Arbol de carpetas listo."

all: $(BIN_A) $(BIN_B) $(BIN_2)
punto1A: $(BIN_A)
punto1B: $(BIN_B)
punto2: $(BIN_2)

# --- Punto 1: el programa escribe en results/ y scripts/ relativos a la
#     carpeta del punto, por eso se ejecuta desde ahi.
run1A: $(BIN_A)
	cd $(PA) && ./bin/punto1A
sim1A: $(BIN_A)
	cd $(PA) && printf '0\n' | ./bin/punto1A
conv1A: $(BIN_A)
	cd $(PA) && printf '1\n' | ./bin/punto1A

run1B: $(BIN_B)
	cd $(PB) && ./bin/punto1B
sim1B: $(BIN_B)
	cd $(PB) && printf '0\n' | ./bin/punto1B
conv1B: $(BIN_B)
	cd $(PB) && printf '1\n' | ./bin/punto1B

# --- Punto 2: el codigo usa rutas "../results/" y "../scripts/", por eso se
#     ejecuta desde punto2_motor_molecular/bin y se crean las carpetas antes.
define RUN2_PREP
mkdir -p $(P2)/results/animations $(P2)/scripts
endef

run2: $(BIN_2)
	$(RUN2_PREP)
	cd $(P2)/bin && ./motor
demo2-det: $(BIN_2)
	$(RUN2_PREP)
	cd $(P2)/bin && printf '%s\n' $(DEMO_DET) | ./motor
demo2-est: $(BIN_2)
	$(RUN2_PREP)
	cd $(P2)/bin && printf '%s\n' $(DEMO_EST) | ./motor

informe:
	cd $(INF) && pdflatex -interaction=nonstopmode examen2_informe.tex
	-cd $(INF) && bibtex examen2_informe
	cd $(INF) && pdflatex -interaction=nonstopmode examen2_informe.tex
	cd $(INF) && pdflatex -interaction=nonstopmode examen2_informe.tex

clean:
	rm -rf $(PA)/build $(PA)/bin $(PB)/build $(PB)/bin $(P2)/build $(P2)/bin
	@echo "build/ y bin/ eliminados (results/ y scripts/ no se tocan)."

clean-results:
	for p in $(PA) $(PB) $(P2); do \
		find $$p/results -type f \( -name '*.dat' -o -name '*.png' -o -name '*.gif' \) -delete 2>/dev/null; \
	done
	@echo "Resultados generados eliminados (scripts/ no se toca)."
