/**
 * @author   Camilo Huertas, Isabel Nieto
 * @date     2025-04-27
 * @version  1.0.0
 * @license  MIT
 */

#include "../include/svlaplaceEquation.h"
#include <iostream>
#include <vector>

using namespace std;


int main() {
	double base, altura, V_izq, V_base, V_escalera, error, lambda, metodo;
	int nx, ny;

	solicitarDatos(base, altura, nx, ny, V_izq, V_base, V_escalera, error, lambda, metodo);

	double dx = base / (nx - 1);
	double dy = altura / (ny - 1);

	vector<vector<double>> matriz;
	inicializarMatriz(matriz, nx, ny, V_izq, V_base, V_escalera);

	if (metodo == 1){
	coeficiente_trapecio(matriz, nx, ny, error, lambda, V_escalera);

	} else if (metodo ==2) {
	coeficiente_simpson(matriz, nx, ny, error, lambda, V_escalera);
	}

	calcularpotencial();
	guardarDatos_eigen(matriz, nx, ny, dx, dy);
	graficarDatos();

	return 0;
}

