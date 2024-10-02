#ifndef CHOLESKY_H
#define CHOLESKY_H
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
using namespace std;
void printMatrix(vector<vector<double>> matrix, int n);
vector<vector<double>> cholesky(vector<vector<double>> matrix, int n);
vector<vector<double>> mat_mul(vector<vector<double>> A, 
vector<vector<double>> B, int n);
double norm(vector<vector<double>> A, vector<vector<double>> B, int n);


#endif
