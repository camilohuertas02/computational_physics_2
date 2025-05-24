/**
 * @file     svlaplaceEquation.cpp
 * @brief    Implementación de funciones para resolver numericamente la ecuacion de laplace 
 * @author   Camilo Huertas, Isabel Nieto
 * @date     2025-05-8
 * @version  1.0.0
 * @license  MIT
 */

#include "../include/svlaplaceEquation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip> // ✅ para setw y setprecision
#include <limits> // ⬅️ Agregar esta línea si no está

using namespace std;

void solicitarDatos(int &Nfourier, int &nx, int &ny, int &metodo, double &radio) {
	cout << "\n |==============================================| \n";
	cout << "\n |========= CONFIGURACIÓN DEL PROBLEMA =========| \n";
	cout << "\n |==============================================| \n";

	cout << "Ingrese el radio del circulo: " << endl;
	cin >> radio;

	cout << "Ingrese metodo para integrar 1. simpson, 2. trapecio: "<< endl;
	cin >> metodo;

	cout << "Ingrese el número de nodos en x (entero >= 3): "<< endl;
	cin >> nx;

	cout << "Ingrese el número de nodos en y (entero >= 3): "<< endl;
	cin >> ny;

	cout << "Ingrese cuantos N terminos de Fourier quiere: " << endl;
	cin >> Nfourier;

	cout << "\n✅ Todos los datos fueron ingresados correctamente.\n";
}

void inicializarMatriz(vector<vector<double>> &matriz, int nx, int ny, int radio, double potencial_0) {
	matriz.resize(nx, vector<double>(ny, 0.0));

}

// Constante π portátil
static constexpr double PI = acos(-1.0);

double compute_Cn_simpson(int n, double a, double V0, int M_simp) {
	const double b = PI/2.0;
	double h = b / M_simp;

	// Integrando I_n = ∫₀^{π/2} V0·sin(2θ)·sin((2n-1)θ) dθ
	auto f = [&](double theta) {
		return V0 * std::sin(2.0 * theta)
			* std::sin((2 * n - 1) * theta);
	};

	// Simpson: extremos
	double s = f(0.0) + f(b);
	// Simpson: internos
	for (int j = 1; j < M_simp; ++j) {
		double factor = (j % 2 == 0 ? 2.0 : 4.0);
		s += factor * f(j * h);
	}

	double I_n = s * h / 3.0;
	// C_n = (4 / (π * a^{2n-1})) · I_n
	return (4.0 / (PI * std::pow(a, 2 * n - 1))) * I_n;
}

void coeficiente_simpson(int nx, int ny, double radio, int Nfourier, vector<double> &coeficientes, double potencial_0) {
	int M_simp = 200;    // elige tú la precisión (debe ser par)

	for (int n = 1; n <= Nfourier; ++n) {
		coeficientes[n] = compute_Cn_simpson(
				n,           // índice del modo
				radio,       // a
				potencial_0, // V0
				M_simp
				);
	}
}


// Calcula C_n usando la regla del trapecio directamente
double compute_Cn_trapezoid(int n, double a, double V0, int M_trap) {
	const double b = PI / 2.0;
	double h = b / M_trap;

	// Integrando I_n = ∫₀^{π/2} V0·sin(2θ)·sin((2n-1)θ) dθ
	auto f = [&](double theta) {
		return V0 * sin(2.0 * theta)
			* sin((2 * n - 1) * theta);
	};

	// Trapecio: suma de extremos
	double s = 0.5 * (f(0.0) + f(b));
	// Trapecio: internos
	for (int j = 1; j < M_trap; ++j) {
		s += f(j * h);
	}

	double I_n = s * h;
	// C_n = (4 / (π * a^{2n-1})) · I_n
	return (4.0 / (PI * pow(a, 2 * n - 1))) * I_n;
}

void coeficiente_trapecio(int nx, int ny, double radio, int Nfourier, vector<double> &coeficientes, double potencial_0) {
	int M_trap = 200;    // elige tú la precisión

	for (int n = 1; n <= Nfourier; ++n) {
		coeficientes[n] = compute_Cn_trapezoid(
				n,            // índice del modo
				radio,        // a
				potencial_0,  // V0
				M_trap
				);
	}
}


void calcularPotencial(vector<vector<double>> &matriz, int nx, int ny, double dx, double dy, double radio, vector<double> coeficientes, int Nfourier){

	for (int i = 0; i < (nx- 1) ; i++){	
		double x = i * dx;
		for (int j = 0; j < (ny - 1); j++){
			double y = j * dy;
			double r = hypot(x, y);	

			if (r > radio){
				matriz[i][j] = 0;
			} else {

				double theta = atan2(y, x);
				double u = 0.0;
				for (int n = 1; n < Nfourier; ++n) {
					int m = 2 * n - 1;
					u += coeficientes[n] * pow(r, m) * sin(m * theta);
				}
				matriz[i][j] = u;
			}
		}
	}
}


void guardarDatos(const vector<vector<double>> &matriz, int nx, int ny, double dx, double dy) {
	// 🧽 Borrar archivos previos para evitar confusión visual
	remove("resultado_gnuplot.png");
	remove("resultado_python.png");

	ofstream file("svlaplace.dat");
	if (!file.is_open()) {
		cerr << "❌ Error al abrir 'svlaplace.dat' para escritura.\n";
		return;
	}

	// 📝 Encabezado (comentario que será ignorado por Gnuplot y Python)
	file << "#    x (m)        y (m)        V(x,y) (V)\n";

	file << fixed << setprecision(6);

	// ⚠️ ¡Orden cambiado! Primero recorre j (y), luego i (x)
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			file << setw(12) << i * dx
				<< setw(12) << j * dy
				<< setw(12) << matriz[i][j] << "\n";
		}
		file << "\n";  // 👈 línea en blanco para que Gnuplot detecte filas
	}

	file.close();
	cout << "✅ Resultados guardados en 'svlaplace.dat'.\n";
}



void graficarDatos() {
	string opcion;
	cout << "\n¿Con qué desea visualizar la gráfica? [gnuplot/python]: ";
	getline(cin, opcion); // sin cin.ignore

	if (opcion == "gnuplot") {
		system("gnuplot -persist scripts/plot_svlaplace.gp");
	} else if (opcion == "python") {
		system("python3 scripts/plot_svlaplace.py");
	} else {
		cout << "❌ Opción no reconocida. Escriba 'gnuplot' o 'python'.\n";
	}
}
