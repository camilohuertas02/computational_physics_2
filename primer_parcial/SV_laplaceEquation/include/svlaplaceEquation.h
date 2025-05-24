/**
 * @file     svlaplaceEquation.h
 * @brief    Declaraciones de funciones para solucionar numericamente la ecuación de Laplace y condiciones de frontera.
 * @author   Camilo Huertas, Isabel Nieto
 * @date     2025-05-08
 * @version  1.0.0
 * @license  MIT
 */

#ifndef SVLAPLACEEQUATION_H
#define SVLAPLACEEQUATION_H

#include <vector>
#include <string>


using namespace std;


//void simpson(vector<vector<double>> &matriz, int nx, int ny, double error, double lambda, double V_escalera);

//void trapecio(Eigen::MatrixXd& matriz, int nx, int ny, double error, double lambda, double V_escalera);

//void calcularpotencial();


void solicitarDatos(int &Nfourier, int &nx, int &ny, int &metodo, double &radio);

void inicializarMatriz(vector<vector<double>> &matriz, int nx, int ny, int radio, double potencial_0);

void coeficiente_simpson(int nx, int ny, double radio, int Nfourier, vector<double> &coeficientes, double potencial_0);

void coeficiente_trapecio(int nx, int ny, double radio, int Nfourier, vector<double> &coeficientes, double potencial_0);

void calcularPotencial(vector<vector<double>> &matriz, int nx, int ny, double dx, double dy, double radio, vector<double> coeficientes, int Nfourier);

void guardarDatos(const vector<vector<double>> &matriz, int nx, int ny, double dx, double dy);

void graficarDatos();


#endif // SVLAPLACEEQUATION_H


