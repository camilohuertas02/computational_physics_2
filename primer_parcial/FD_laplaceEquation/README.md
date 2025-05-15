# Proyecto: Solución Numérica de la Ecuación de Laplace

Este programa resuelve la ecuación de Laplace en 2D usando el método de diferencias finitas con Gauss-Seidel.

## Características
- Entrada interactiva de parámetros físicos y térmicos.
- Sobrerrelajación opcional (λ).
- Exportación de resultados a archivo `.dat`.
- Visualización en 3D usando Gnuplot o Python.

## Requisitos
- C++11 o superior
- Gnuplot o Python 3 con matplotlib/numpy

## Compilación y ejecución

```bash
make
make run
```

## Visualización
Se le preguntará si desea graficar con `gnuplot` o `python`. El archivo de datos es `laplace.dat`.

## Archivos

- `include/laplaceEquation.h`: prototipos
- `src/laplaceEquation.cpp`: funciones
- `src/laplaceEquationMain.cpp`: función `main`
- `scripts/plot_laplace.py`: visualización Python
