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

void solicitarDatos(double &base, double &altura, int &nx, int &ny,
		double &V_izq, double &V_base, double &V_escalera,
		double &error, double &lambda, double &metodo);

void validarDatos(double &dato, const string &mensaje, double minimo = 0.0, double maximo = 1e9);

void inicializarMatriz(vector<vector<double>> &matriz, int nx, int ny, double V_izq, double V_base, double V_escalera);

void simpson(vector<vector<double>> &matriz, int nx, int ny, double error, double lambda, double V_escalera);

void trapecio(Eigen::MatrixXd& matriz, int nx, int ny, double error, double lambda, double V_escalera);

void calcularpotencial();

void guardarDatos(const vector<vector<double>> &matriz, int nx, int ny, double dx, double dy);

void graficarDatos();

#endif // SVLAPLACEEQUATION_H


