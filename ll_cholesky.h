#ifndef LL_CHOLESKY_H
#define LL_CHOLESKY_H
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include "my_exception.h"
#include "cholesky.h"
using namespace std;
void printMatrix_rect(vector<vector<double>> A, int m, int n);
vector<vector<double>> mat_mul_rect(vector<vector<double>> A, int m, int n,
		vector<vector<double>> B, int p);
vector<vector<double>> mat_add(vector<vector<double>>& A, int p_0, int p_1,
		int q_0, int q_1, vector<vector<double>> B, 
		int m, int n, int add_sub=1, int overwrite=0);
vector<vector<double>> build_submatrix(vector<vector<double>> A, int m, int
		n, int p_0, int p_1, int q_0, int q_1);
vector<vector<double>> transpose(vector<vector<double>> A, int m, int n);
vector<vector<double>> ll_cholesky(vector<vector<double>> matrix, int n);

#endif
