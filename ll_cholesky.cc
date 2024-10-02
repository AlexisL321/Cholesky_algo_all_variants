#include "ll_cholesky.h"
using namespace std;

//function for printing out rectangular matrices
void printMatrix_rect(vector<vector<double>> A, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << A[i][j] << "  ";
		}
		cout<<endl;
	}
}

//function for rectangular matrix multiplication
vector<vector<double>> mat_mul_rect(vector<vector<double>> A, int m, int n,
		vector<vector<double>> B, int p) {
	//matrix A is of size mxn, matrix B is of size nxp
	vector<vector<double>> result(m, vector<double>(p));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			double sum = 0;
			for (int k = 0; k < n; k++) {
				sum += A[i][k]*B[k][j];
			}
			result[i][j] = sum;
		}
	}
	return result;
}

//function for adding or subtracting two matrices
//or updating part of matrix A with B
//A[p_0:p_1][q_0:q_1] will be updated by B of size mxn
vector<vector<double>> mat_add(vector<vector<double>>& A, int p_0, int
		p_1, int q_0, int q_1, vector<vector<double>> B, int m,
		int n, int add_sub, int overwrite) {
	if (p_1-p_0+1 != m || q_1-q_0+1 != n) {
		throw My_exception("Error: check your argument for dimension "
				"of matrices.");
	}
	vector<vector<double>> result;
	if (overwrite == 0) {
		//not overwriting matrix A
		result = vector<vector<double>>(m, vector<double>(n));
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (add_sub == 1) {
				//matrix addition
				(overwrite == 0) ? result[i][j] = A[i+p_0][j+q_0] + B[i][j] : 
					A[i+p_0][j+q_0] += B[i][j];
			}
			else {
				//matrix subtraction
				(overwrite == 0) ? result[i][j] = A[i+p_0][j+q_0] - B[i][j] : 
					A[i+p_0][j+q_0] -= B[i][j];
			}
		}
	}
	return (overwrite == 0) ? result : A;
	//(overwrite == 0) ? return result : return A;
}

//function to build a submatrix (p_1, q_1 inclusive)
//p is the row and q is the column
vector<vector<double>> build_submatrix(vector<vector<double>> A, int m,
		int n, int p_0, int p_1, int q_0, int q_1) {
	if (p_0 < 0 || p_1 >= m || q_0 < 0 || q_1 >=n) {
		throw My_exception("submatrix length/width is beyond "
				"original matrix");
	}
	vector<vector<double>> result(p_1-p_0+1, vector<double>(q_1-q_0+1));
	for (int i = 0; i < p_1-p_0+1; i++) {
		for (int j = 0; j < q_1-q_0+1; j++) {
			result[i][j] = A[i+p_0][j+q_0];
		}
	}
	return result;
}

//function to transpose matrix
vector<vector<double>> transpose(vector<vector<double>> A, int m, int n){
	vector<vector<double>> result(n, vector<double>(m));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			result[j][i] = A[i][j];
		}
	}
	return result;
}

//left-looking cholesky factorization
vector<vector<double>> ll_cholesky(vector<vector<double>> A, int n) {
	vector<vector<double>> result(n, vector<double>(n));

	for (int j = 0; j < n; j++) {
		if (j != 0) {
			//build submatrix A[j:n,1:j-1] (<-index start from 1)
			vector<vector<double>> sub_A;
			sub_A = build_submatrix(A, n, n, j, (n-1), 0, (j-1));

			//build submatrix A[j, 1:j-1] 9<-index starts from 1)
			vector<vector<double>> vec_A;
			vec_A = build_submatrix(A, n, n, j, j, 0, (j-1));
			vector<vector<double>> vec_A_T;
			vec_A_T = transpose(vec_A, 1, j);

			//multiply sub_A and vec_A_T
			vector<vector<double>> update;
			update = mat_mul_rect(sub_A, n-j, j, vec_A_T, 1);

			/*
			   if (j == 8) {
			   printMatrix_rect(update, n-j, 1);
			   cout<<A[j][j]<<endl;
			   }
			 */

			//update A[j:n,j] with the matrix "update"
			A = mat_add(A, j, n-1, j, j, update, n-j, 1, 0, 1);
		}
		//cout<<"j: "<<j<<endl;
		//cout << A[j][j] <<endl;


		//update diagonal element of A
		A[j][j] = sqrt(A[j][j]);
		//update entries below and above the diagonal
		for (int i = j+1; i < n; i++){
			A[i][j] = A[i][j]/A[j][j];
		}
		for (int i = 0; i < j; i++) {
			A[i][j] = 0;
		}

	}
	return A;
}

/*
int main () {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> rand(0.0,1.0);

	int size = 500;
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

	vector<vector<double>> L = ll_cholesky(A, size);
	vector<vector<double>> LL = mat_mul_rect(L,size, size,
			transpose(L,size,size),size);
	//printMatrix_rect(L, size, size);
	cout<<"norm: "<<norm(A, LL, size) <<endl;
}
*/
