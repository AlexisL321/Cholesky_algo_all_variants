#ifndef LL_BLC_CHO_H
#define LL_BLC_CHO_H

#include "ll_cholesky.h"
#include "block_cho.h"

using namespace std;
vector<vector<double>> update_mat(vector<vector<double>>& A, int p1, int p2,
        int q1, int q2, vector<vector<double>> B, bool is_constant=false,
        int scalar=0);
vector<vector<double>> ll_blc(vector<vector<double>> A, int m, int n,
        int b);

#endif
