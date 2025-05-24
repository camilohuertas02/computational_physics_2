# =============================================
# ============ Makefile: Laplace ==============
# =============================================

# Compilador y opciones
CXX = g++
CXXFLAGS = -Wall -std=c++11

# Directorios
INCLUDE_DIR = include
SRC_DIR     = src
BIN         = laplaceEquation

# Archivos fuente
SOURCES = $(SRC_DIR)/laplaceEquation.cpp $(SRC_DIR)/laplaceEquationMain.cpp

# Regla por defecto
all: $(BIN)

# Cómo compilar el ejecutable
$(BIN): $(SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -o $(BIN) $(SOURCES)

# Ejecutar el programa
run: $(BIN)
	./$(BIN)

# Limpiar archivos generados
clean:
	rm -f $(BIN) *.o laplace.dat resultado_*.png resultado_*.gp

# Regla completa
full: clean all run

.PHONY: all clean run full


