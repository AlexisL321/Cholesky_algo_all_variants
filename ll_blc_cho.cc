#include "ll_blc_cho.h"
using namespace std;

// ll-cholesky block factorization 
//matrix size m, block size n, number of block b(b=4 means there are 4
//blocks 0, 1, 2, 3
vector<vector<double>> ll_blc(vector<vector<double>> A, int m, int n, 
		int b) {
	vector<vector<double>> result(n, vector<double>(n));
	for (int j = 0; j < b; j++) {
		if (j != 0) {
			//build submatrix A[j*n+1:b*n,1:j*n] (<-index start from 1)
			vector<vector<double>> sub_A;
			sub_A = build_submatrix(A, m, m, j*n, (b*n-1), 0, (j*n-1));
			//cout<<"j:"<<j<<endl;//TODO
			//printMatrix_rect(sub_A, (b-j)*n, j*n);

			//build submatrix A[j*n+1, 1:j*n] (<-index starts from 1)
			vector<vector<double>> vec_A;
			vec_A = build_submatrix(A, m, m, j*n, (j+1)*n-1, 0, (j*n-1));
			vector<vector<double>> vec_A_T;
			vec_A_T = transpose(vec_A, n, j*n); 

			//multiply sub_A and vec_A_T
			vector<vector<double>> update;
			//TODO: can parallelize the matrix multiplication
			update = mat_mul_rect(sub_A, (b-j)*n, j*n, vec_A_T, n); 

			//update A[j*n:b*n-1, j*n:(j+1)*n-1] with the matrix "update"
			A = mat_add(A, j*n, b*n-1, j*n, (j+1)*n-1, update, 
					(b-j)*n, n, 0, 1); 
		}   
		//cout<<"j: "<<j<<endl;
		//cout << A[j][j] <<endl;


		//update diagonal block of A
		vector<vector<double>> diag_A = build_submatrix(A, m, m, 
				j*n, (j+1)*n-1, j*n, (j+1)*n-1);
		vector<vector<double>> diag_L = ll_cholesky(diag_A, n);
		//update matrix A diagonal block
		update_mat(A, j*n, (j+1)*n-1, j*n, (j+1)*n-1, diag_L);
		//update blocks below and above the diagonal
		//first get the inverse of diag_L
		vector<vector<double>> diag_L_inv = lower_inverse(diag_L, n);

		//updating blocks below the diagonal block
		for (int i = j+1; i < b; i++) {
			vector<vector<double>> copy = build_submatrix(A, m, m, i*n, 
					(i+1)*n-1, j*n, (j+1)*n-1);
			//printMatrix_rect(copy, n, n);
			//cout<<i<<" "<<j<<endl;//TODO
			//transpose of diag_L_inv instead of just dia_L_inv
			copy = mat_mul_rect(copy, n, n, transpose(diag_L_inv, n, n), n);
			update_mat(A, i*n, (i+1)*n-1, j*n, (j+1)*n-1, copy);

		}

		//updating blocks above diagonal block
		for (int i = 0; i < j; i++) {
			//placeholder matrix
			vector<vector<double>> ph(n, vector<double>(n));
			update_mat(A, i*n, (i+1)*n-1, j*n, (j+1)*n-1, ph, true, 0);
		}

	}
	return A;

}

//function to update a part of a matrix with either another matrix B
//or with a constant scalar, is_constant is defaulted to be false, scalar
//is defaulted to be 0. the size of B must equal to the size of the
//submatrix of A to be updated. p1, p2 --> rows, q1, q2 --> columns,
//inclusive of p1, p2, q1, q2
vector<vector<double>> update_mat(vector<vector<double>>& A, int p1, int p2,
		int q1, int q2, vector<vector<double>> B, bool is_constant, 
		int scalar) {
	//dimension of B
	for (int i = p1; i <= p2; i++) {
		for (int j = q1; j <= q2; j++) {
			A[i][j] = (is_constant) ? scalar : B[i-p1][j-q1];
		}
	}
	return A;
}

int main() {
	/*
	//test1 test the update_mat function -- passed
	vector<vector<double>> test_A = {{1,2,3},{4,5,6},{7,8,9}};
	vector<vector<double>> test_B = {{0,0},{0,0}};
	update_mat(test_A, 0, 1, 0, 1, test_B, true, 1);
	printMatrix_rect(test_A, 3, 3);
	 */
	/*
	//test with a simple matrix -- passed
	vector<vector<double>> test_L = {{1,0,0,0},{2,5,0,0},{3,6,8,0},{4,7,9,10}};
	vector<vector<double>> test_A = {{1,2,3,4},{2,29,36,43},{3,36,109,126},
	{4,43,126,246}};
	vector<vector<double>> result_L = ll_blc(test_A, 4, 1, 4);
	printMatrix_rect(result_L,4, 4);
	 */
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> rand(0.0,1.0);

	int size = 500;
	int block_size = 10;
	int b = floor(size/block_size);
	vector<vector<double>> A(size, vector<double>(size));
	for (int i = 0; i < size; i++) {
		for (int j =0; j < size; j++) {
			A[i][j] = rand(gen);
		}
	}
	//make A symmetric
	A = mat_add(A,0,size-1,0,size-1,transpose(A,size,size),size,size,1,1);
	//cout<<"original matrix A: "<<endl;
	//printMatrix_rect(A,size,size);

	//make A positive definite
	for (int i = 0; i < size; i++) {
		A[i][i]+=size;
	}
	//test
	/**
	vector<vector<double>> L_tst =
	{{1,0,0,0,0,0},{1,2,0,0,0,0},{1,2,3,0,0,0},{1,2,3,4,0,0},{1,2,3,4,5,0},
	{1,2,3,4,5,6}};
	A = mat_mul_rect(L_tst, 6, 6,
	transpose(L_tst, 6, 6), 6);
	**/

	vector<vector<double>> L = ll_blc(A, size, block_size, b);

	vector<vector<double>> LL = mat_mul_rect(L,size, size,
			transpose(L,size,size),size);
	printMatrix_rect(L, size, size);
	cout<<"norm: "<<norm(A, LL, size) <<endl;

}
