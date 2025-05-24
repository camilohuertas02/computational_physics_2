/**
 * @file     laplaceEquation.cpp
 * @brief    Implementación de funciones para resolver la ecuación de Laplace usando el método de Gauss-Seidel y visualizar resultados con Gnuplot o Python.
 * @author   John H. Díaz Forero
 * @date     2025-04-27
 * @version  1.0.0
 * @license  MIT
 */

#include "../include/laplaceEquation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip> // ✅ para setw y setprecision
#include <limits> // ⬅️ Agregar esta línea si no está


using namespace std;


void validarDatos(double &dato, const string &mensaje, double minimo, double maximo) {
    string entrada;
    stringstream ss;

    do {
        getline(cin, entrada);
        ss.clear();
        ss.str(entrada);

        if (!(ss >> dato) || dato < minimo || dato > maximo) {
            cout << "❌ Error: " << mensaje << " (rango permitido: " << minimo << " a " << maximo << "). Intente nuevamente: ";
        } else {
            break;
        }

    } while (true);
}

void dibujarRectangulo(double base, double altura) {
    const int MAX_ANCHO = 40;  // caracteres (base)
    const int MAX_ALTO  = 10;  // líneas (altura)
  

    int ancho = min(MAX_ANCHO, static_cast<int>(round(base * 10)));
    int alto  = min(MAX_ALTO,  static_cast<int>(round(altura * 10)));

    cout << "\n|==============================================| \n";
    cout << "|      REPRESENTACIÓN DEL DOMINIO RECTANGULAR  |\n";
    cout << "|==============================================| \n\n";

   // Línea superior
    cout << "  +";
    for (int i = 0; i < ancho; ++i) cout << "-";
    cout << "+\n";

    // Lados y texto de altura en la mitad
    for (int j = 0; j < alto; ++j) {
        cout << "  |";
        for (int i = 0; i < ancho; ++i) cout << " ";
        cout << "|";

        if (j == alto / 2)
            cout << "   ← altura (" << altura << " m)";
        cout << "\n";
    }

    // Línea inferior
    cout << "  +";
    for (int i = 0; i < ancho; ++i) cout << "-";
    cout << "+\n";

    // Flecha base
    int arrow_pos = ancho / 2;
    cout << string(2 + arrow_pos, ' ') << "↑\n";
    cout << string(2 + arrow_pos - 4, ' ') << "base (" << base << " m)\n\n";
}
void solicitarDatos(double &base, double &altura, int &nx, int &ny,
                    double &T_izq, double &T_der, double &T_sup, double &T_inf,
                    double &error, double &lambda) {
    double entradaTemporal;
  cout << "\n |==============================================| \n";
  cout << "\n |========= CONFIGURACIÓN DEL PROBLEMA =========| \n";
  cout << "\n |==============================================| \n";

    cout << "Ingrese la base del rectángulo (m > 0): ";
    validarDatos(base, "La base debe ser mayor que 0", 0.01);

    cout << "Ingrese la altura del rectángulo (m > 0): ";
    validarDatos(altura, "La altura debe ser mayor que 0", 0.01);
    dibujarRectangulo(base, altura); 
    cout << "Ingrese el número de nodos en x (entero >= 3): ";
    validarDatos(entradaTemporal, "Número de nodos en x inválido", 3);
    nx = static_cast<int>(entradaTemporal);

    cout << "Ingrese el número de nodos en y (entero >= 3): ";
    validarDatos(entradaTemporal, "Número de nodos en y inválido", 3);
    ny = static_cast<int>(entradaTemporal);

    cout << "Temperatura del borde izquierdo (°C): ";
    validarDatos(T_izq, "Temperatura inválida");

    cout << "Temperatura del borde derecho (°C): ";
    validarDatos(T_der, "Temperatura inválida");

    cout << "Temperatura del borde superior (°C): ";
    validarDatos(T_sup, "Temperatura inválida");

    cout << "Temperatura del borde inferior (°C): ";
    validarDatos(T_inf, "Temperatura inválida");

    cout << "Ingrese el error permitido para la convergencia (%): ";
    validarDatos(error, "El error debe estar entre 0 y 100", 0.01, 100.0);

    cout << "Ingrese λ para sobrerrelajación (1 ≤ λ ≤ 2) (Recomendación: si no está seguro, ingrese λ = 1.0 (sin   sobrerrelajación)): ";
    
    validarDatos(lambda, "λ debe estar entre 1.0 y 2.0", 1.0, 2.0);
	cin.ignore(numeric_limits<streamsize>::max(), '\n');  // 👉 Agrega esto aquí
    cout << "\n✅ Todos los datos fueron ingresados correctamente.\n";
}

void inicializarMatriz(vector<vector<double>> &matriz, int nx, int ny,
                       double T_izq, double T_der, double T_sup, double T_inf) {
    matriz.resize(nx, vector<double>(ny, 0.0));

    for (int i = 0; i < nx; i++) {
        matriz[i][0] = T_inf;
        matriz[i][ny - 1] = T_sup;
    }

    for (int j = 0; j < ny; j++) {
        matriz[0][j] = T_izq;
        matriz[nx - 1][j] = T_der;
    }
}

void gaussSeidel(vector<vector<double>> &matriz, int nx, int ny, double error, double lambda) {
    double maxError;
    int iter = 0;

    do {
        maxError = 0.0;

        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                double T_old = matriz[i][j];
                double T_new = 0.25 * (matriz[i+1][j] + matriz[i-1][j] +
                                       matriz[i][j+1] + matriz[i][j-1]);
                matriz[i][j] = lambda * T_new + (1 - lambda) * T_old;

                double diff = fabs((matriz[i][j] - T_old) / matriz[i][j]) * 100.0;
                if (diff > maxError) maxError = diff;
            }
        }

        iter++;

    } while (maxError > error);

    cout << "✅ Convergió en " << iter << " iteraciones.\n";
}

void guardarDatos(const vector<vector<double>> &matriz, int nx, int ny, double dx, double dy) {
    // 🧽 Borrar archivos previos para evitar confusión visual
    remove("resultado_gnuplot.png");
    remove("resultado_python.png");

    ofstream file("laplace.dat");
    if (!file.is_open()) {
        cerr << "❌ Error al abrir 'laplace.dat' para escritura.\n";
        return;
    }

    // 📝 Encabezado (comentario que será ignorado por Gnuplot y Python)
    file << "#    x (m)        y (m)        T(x,y) (°C)\n";

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
    cout << "✅ Resultados guardados en 'laplace.dat'.\n";
}
void graficarDatos() {
    string opcion;
    cout << "\n¿Con qué desea visualizar la gráfica? [gnuplot/python]: ";
    getline(cin, opcion); // sin cin.ignore

    if (opcion == "gnuplot") {
        system("gnuplot -persist scripts/plot_laplace.gp");
    } else if (opcion == "python") {
        system("python3 scripts/plot_laplace.py");
    } else {
        cout << "❌ Opción no reconocida. Escriba 'gnuplot' o 'python'.\n";
    }
}
