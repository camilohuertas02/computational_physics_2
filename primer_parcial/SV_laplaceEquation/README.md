# svLaplaceEquation

Implementación en C++ para la resolución numérica de la ecuación de Laplace en coordenadas polares, usando desarrollos en series de Fourier y métodos de integración (Simpson o Trapecio).

---

## Índice

1. [Descripción](#descripción)  
2. [Características](#características)  
3. [Requisitos](#requisitos)  
4. [Estructura del proyecto](#estructura-del-proyecto)  
5. [Compilación](#compilación)  
6. [Uso](#uso)  
   - [Ejecutar el programa](#ejecutar-el-programa)  
   - [Opciones de visualización](#opciones-de-visualización)  
7. [Detalles de los parámetros](#detalles-de-los-parámetros)  
8. [Formato de salida](#formato-de-salida)  
9. [Ejemplo](#ejemplo)  
10. [Autoría y licencia](#autoría-y-licencia)  

---

## Descripción

Este proyecto resuelve la ecuación de Laplace dentro de un dominio circular de radio `a`, imponiendo un potencial de contorno `V0` y desplegando la solución mediante un desarrollo en series de Fourier:

\[
V(r,\theta)=\sum_{n=1}^{N}C_n\,r^{2n-1}\,\sin((2n-1)\theta)
\]

Los coeficientes \(C_n\) se calculan numéricamente mediante la regla de Simpson o del Trapecio aplicadas a la integral:

\[
C_n = \frac{4}{\pi\,a^{2n-1}} \int_0^{\pi/2} V_0\,\sin(2\theta)\,\sin((2n-1)\theta)\,d\theta.
\]

---

## Características

- Cálculo de coeficientes \(C_n\) con dos métodos de integración: Simpson y Trapecio.
- Armado de la malla rectangular en coordenadas cartesianas cubriendo el interior del círculo.
- Exportación de resultados a fichero de datos `svlaplace.dat` compatible con Gnuplot o Python.
- Opción para graficar directamente los resultados vía Gnuplot o script Python.

---

## Requisitos

- **Compilador C++** con soporte C++11 o superior (g++, clang++).
- Librerías estándar de C++ (`<iostream>`, `<vector>`, `<cmath>`, `<fstream>`, `<iomanip>`, `<limits>`).
- **Gnuplot** (si se desea visualización via Gnuplot).
- **Python 3** con `matplotlib` y `numpy` (si se desea visualización via script Python).

---

## Estructura del proyecto

```plain
project-root/
├── include/
│   └── svlaplaceEquation.h        # Declaraciones de funciones
├── src/
│   ├── svlaplaceEquation.cpp      # Implementación de funciones
│   └── svlaplaceEquationMain.cpp  # Función main y flujo principal
├── scripts/
│   ├── plot_svlaplace.gp          # Script Gnuplot
│   └── plot_svlaplace.py          # Script Python
└── README.md                      # Este archivo
```

---

## Compilación

Desde la raíz del proyecto, por ejemplo usando `g++`:

```bash
mkdir -p build && cd build
cmake ..                   # si usas CMakeLists
make                       # o:
g++ -std=c++11 -I../include ../src/svlaplaceEquation.cpp ../src/svlaplaceEquationMain.cpp -o svlaplace
```

O compilación directa sin CMake:

```bash
g++ -std=c++11 -Iinclude src/svlaplaceEquation.cpp src/svlaplaceEquationMain.cpp -o svlaplace
```

---

## Uso

### Ejecutar el programa

```bash
./svlaplace
```

El programa pedirá interactivamente:

1. **Radio del círculo** (`double`).  
2. **Método de integración**: `1` para Simpson, `2` para Trapecio.  
3. **Número de nodos en X** (`int >= 3`).  
4. **Número de nodos en Y** (`int >= 3`).  
5. **Número de términos de Fourier** (`int > 0`).  

Tras introducirlos, se calculará la distribución de potencial y se guardarán los datos en `svlaplace.dat`.

### Opciones de visualización

Después, se ofrece la opción de graficar con:

- **gnuplot**: ejecuta el script `scripts/plot_svlaplace.gp`.  
- **python**: ejecuta `scripts/plot_svlaplace.py`.  

---

## Detalles de los parámetros

- `Nfourier`: número de modos de Fourier a incluir (modos impares).  
- `nx`, `ny`: discretización en X e Y.  
- `radio`: radio del dominio circular.  
- `potencial_0`: valor de potencial en la frontera (V0).  

---

## Formato de salida

Archivo `svlaplace.dat` con columnas:

```
#    x (m)        y (m)        V(x,y) (V)
x_0    y_0    V(x_0,y_0)
...
```

Separado por espacios y bloques en blanco para facilitar Gnuplot.

---

## Ejemplo

```bash
$ ./svlaplace
Ingrese el radio del circulo:
> 1.0
Ingrese metodo para integrar 1. simpson, 2. trapecio:
> 1
Ingrese el número de nodos en x (entero >= 3):
> 100
Ingrese el número de nodos en y (entero >= 3):
> 100
Ingrese cuantos N terminos de Fourier quiere:
> 10

✅ Resultados guardados en 'svlaplace.dat'.

¿Con qué desea visualizar la gráfica? [gnuplot/python]:
> python
```

---

## Autoría y licencia

- **Autores**: Camilo Huertas, Isabel Nieto  
- **Fecha**: 8 de mayo de 2025  
- **Versión**: 1.0.0  
- **Licencia**: MIT. Consulte el fichero `LICENSE` para más detalles.  
