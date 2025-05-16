/**
 * @file     laplaceEquation.h
 * @brief    Declaraciones de funciones para resolver la ecuación de Laplace con Gauss-Seidel y condiciones de frontera.
 * @author   Camilo Huertas, Isabel Nieto
 * @date     2025-05-08
 * @version  1.0.0
 * @license  MIT
 */

#ifndef LAPLACEEQUATION_H
#define LAPLACEEQUATION_H

#include <vector>
#include <string>
#include <Eigen/Dense>

using namespace std;

void solicitarDatos(double &base, double &altura, int &nx, int &ny,
		double &V_izq, double &V_base, double &V_escalera,
		double &error, double &lambda, double &metodo);

void validarDatos(double &dato, const string &mensaje, double minimo = 0.0, double maximo = 1e9);

void inicializarMatriz(vector<vector<double>> &matriz, int nx, int ny, double V_izq, double V_base, double V_escalera);

void inicializarMatriz_eigen(Eigen::MatrixXd& matriz , int nx, int ny, double V_izq, double V_base, double V_escalera);

void gaussSeidel(vector<vector<double>> &matriz, int nx, int ny, double error, double lambda, double V_escalera);

void gaussSeidel_eigen(Eigen::MatrixXd& matriz, int nx, int ny, double error, double lambda, double V_escalera);

void guardarDatos(const vector<vector<double>> &matriz, int nx, int ny, double dx, double dy);

void guardarDatos_eigen(const Eigen::MatrixXd& matriz, int nx, int ny, double dx, double dy);

void graficarDatos();

#endif // LAPLACEEQUATION_H

