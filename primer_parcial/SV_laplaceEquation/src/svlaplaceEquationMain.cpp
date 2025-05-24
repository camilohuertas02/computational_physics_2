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
	double radio;
	double potencial_0 = 50; 
	int nx, ny, Nfourier, metodo;

	solicitarDatos(Nfourier, nx, ny, metodo, radio);

	double dx = radio / (nx - 1);
	double dy = radio / (ny - 1);

	vector<vector<double>> matriz;
	inicializarMatriz(matriz, nx, ny, radio, potencial_0);

	vector<double> coeficientes(Nfourier + 1);

	if (metodo == 1){
		coeficiente_simpson(nx, ny, radio, Nfourier, coeficientes, potencial_0);

	} else if (metodo ==2) {
		coeficiente_trapecio(nx, ny, radio, Nfourier, coeficientes, potencial_0);	
	}

	calcularPotencial(matriz, nx, ny, dx, dy, radio, coeficientes, Nfourier);

	guardarDatos(matriz, nx, ny, dx, dy);

	graficarDatos();

	return 0;
}

