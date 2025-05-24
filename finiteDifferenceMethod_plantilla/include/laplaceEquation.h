/**
 * @file     laplaceEquation.h
 * @brief    Declaraciones de funciones para resolver la ecuación de Laplace con Gauss-Seidel y condiciones de frontera.
 * @author   
 * @date     2025-04-27
 * @version  1.0.0
 * @license  MIT
 */

#ifndef LAPLACEEQUATION_H
#define LAPLACEEQUATION_H

#include <vector>
#include <string>

using namespace std;

void solicitarDatos(double &base, double &altura, int &nx, int &ny,
                    double &T_izq, double &T_der, double &T_sup, double &T_inf,
                    double &error, double &lambda);

void validarDatos(double &dato, const string &mensaje, double minimo = 0.0, double maximo = 1e9);

void inicializarMatriz(vector<vector<double>> &matriz, int nx, int ny,
                        double T_izq, double T_der, double T_sup, double T_inf);

void gaussSeidel(vector<vector<double>> &matriz, int nx, int ny, double error, double lambda);

void guardarDatos(const vector<vector<double>> &matriz, int nx, int ny, double dx, double dy);

void graficarDatos();

#endif // LAPLACEEQUATION_H

