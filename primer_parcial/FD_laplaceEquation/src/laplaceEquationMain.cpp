/**
 * @author
 * @date     2025-04-27
 * @version  1.0.0
 * @license  MIT
 */

#include "../include/laplaceEquation.h"
#include <iostream>
#include <vector>

using namespace std;


int main() {
    double base, altura, T_izq, T_der, T_sup, T_inf, error, lambda;
    int nx, ny;

    solicitarDatos(base, altura, nx, ny, T_izq, T_der, T_sup, T_inf, error, lambda);

    double dx = base / (nx - 1);
    double dy = altura / (ny - 1);

    vector<vector<double>> matriz;
    inicializarMatriz(matriz, nx, ny, T_izq, T_der, T_sup, T_inf);

    gaussSeidel(matriz, nx, ny, error, lambda);

    guardarDatos(matriz, nx, ny, dx, dy);

    graficarDatos();

    return 0;
}
