//#include <iostream>
//#include <cmath>
//#include <random>
//#include <vector>
#include "cholesky.h"

using namespace std;


void printMatrix(vector<vector<double>> matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

vector<vector<double>> cholesky(vector<vector<double>> matrix, int n) {
	//i is iterating through columns
	for (int i = 0; i < n; i++) {
		//modifying the first column
		for (int j = 0; j < n; j++) {
			if (i == j) {
			matrix[i][j] = sqrt(matrix[i][j]);
			continue;
			}
			if (j < i) {
			matrix[j][i] = 0;
			continue;
			}
			matrix[j][i] = matrix[j][i] / matrix[i][i];
		}

		//modifying the rest of the matrix
		for (int j = i; j < n; j++) {
			//j is iterating through column after i
			if (i == j) continue;
			for (int k = i; k < n; k++) {
				//k is iterating through rows after i
				if (k == i) {
				matrix[k][j] = 0;
				continue;
				}
				matrix[k][j] = matrix[k][j] - matrix[k][i]*matrix[j][i];
			}
		}
	}
return matrix;

}

vector<vector<double>> mat_mul(vector<vector<double>> A, 
		vector<vector<double>> B, int n) {
	vector<vector<double>> result(n, vector<double>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0;
			for (int k = 0; k < n; k++) {
				sum += A[i][k]*B[k][j];
			}
			result[i][j] = sum;
		}
	}
	return result;
}

double norm(vector<vector<double>> A, vector<vector<double>> B, int n) {
	double result = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			result += abs(A[i][j] - B[i][j]);
		}
	}
	return result/(n*n);
}

/*
int main() {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> rand(0.0, 1.0);

	int size = 10000;
	vector<vector<double>> A(size, vector<double>(size));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			A[i][j] = 0.5 * (rand(gen) + rand(gen));
			A[i][j] = A[j][i];
		}
	}

	for (int i = 0; i < size; i++) {
		A[i][i] += 10;
	}

	cout << "A:" << endl;
	printMatrix(A, size);
	vector<vector<double>> alpha(size, vector<double>(size));
	alpha = cholesky(A, size);

	cout << "Cholesky alpha:" << endl;
	printMatrix(alpha, size);
	vector<vector<double>> new_A (size, vector<double>(size));
	vector<vector<double>> alpha_T(size, vector<double>(size));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++){
			alpha_T[i][j] = alpha[j][i];
		}
	}

	new_A = mat_mul(alpha, alpha_T, size);
	cout<<"norm: "<< norm(A, new_A, size)<<endl;
	return 0;
}
*/

