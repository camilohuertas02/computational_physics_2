# ========================================================
# ============ Makefile: LaplaceIntegracion ==============
# ========================================================

# Compilador y opciones
CXX = g++
CXXFLAGS = -Wall -std=c++11 -g

# Directorios
INCLUDE_DIR = include
SRC_DIR     = src
BIN         = svlaplaceEquation

# Archivos fuente
SOURCES = $(SRC_DIR)/svlaplaceEquation.cpp $(SRC_DIR)/svlaplaceEquationMain.cpp

# Regla por defecto
all: $(BIN)

# Cómo compilar el ejecutable
#$(BIN): $(SOURCES)
#	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -o $(BIN) $(SOURCES)

$(BIN): $(SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -o $(BIN) $(SOURCES)

# Ejecutar el programa
run: $(BIN)
	./$(BIN)

# Limpiar archivos generados
clean:
	rm -f $(BIN) *.o svlaplace.dat resultado_*.png resultado_*.gp

# Regla completa
full: clean all run

.PHONY: all clean run full


