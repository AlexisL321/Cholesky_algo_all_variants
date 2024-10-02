#ifndef BLOCK_CHOL_H
#define BLOCK_CHOL_H
#include "cholesky.h"
//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <random>

vector<vector<double>> lower_inverse(vector<vector<double>> A, int n);

vector<vector<double>> block_chol(vector<vector<double>> copy, int n,
int m, int num_block);


#endif
