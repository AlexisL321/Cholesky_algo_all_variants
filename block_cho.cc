//#include "cholesky.h"
//#include <random>
//#include <cmath>
//#include <iostream>
#include "block_cho.h"

//a function to calculate the inverse of a lower triangular matrix
vector<vector<double>> lower_inverse(vector<vector<double>> A, int n) {
	vector<vector<double>> inverse(n, vector<double>(n));
	//copy the original matrix
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			inverse[i][j] = 0;
		}
	}

	//calculate the inverse
	for (int i = 0; i < n; i++) {
		//solve the linear system of equations n times
		vector<double> b(n);
		for (int j = 0; j < n; j++) {
			if (j == i) {
				b[j] = 1;
				continue;
			}
			b[j] = 0;
		}

		for (int j = 0; j < n; j++) {
			//j is the jth row of the matrix
			if (j == 0) {
				inverse[j][i] = b[j]/A[j][j];
				continue;
			}

			double sum = b[j];
			for (int k = 0; k < j; k++){
				//k iterates through each entry in jth row
				sum = sum - A[j][k] * inverse[k][i];
			}
			inverse[j][i] = sum / A[j][j];
		}
	}
	return inverse;
}


vector<vector<double>> block_chol(vector<vector<double>> copy, 
		int n, int m, int num_block) {//m is the block size
									  //num_block is the number of block matrices in a row
	/**
	vector<vector<double>> copy (n, vector<double>(n));

	//copy the input matrix
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			copy[i][j] = matrix[i][j];
		}
	}
	**/

	for (int i = 0; i < num_block; i++) {
		//i iterates through the blocks horizontally
		vector<vector<double>> chol_block (m, vector<double>(m));
		for (int j = 0; j < num_block; j++) {
			//j iterates through the blocks vertically

			//for blocks along the diagonal
			if (i == j) {
				//copy this block matrix
				for (int k = 0; k < m; k++) {
					for (int p = 0; p < m; p++) {
						chol_block[k][p] = copy[i*m+k][j*m+p];
					}
				}
				//do cholesky factorization on this matrix
				chol_block = cholesky(chol_block, m);


				//update the element in original matrix (i.e. copy)
				for (int k = 0; k < m; k++) {
					for (int p = 0; p < m; p++) {
						copy[i*m+k][j*m+p] = chol_block[k][p];
					}
				}
				continue;
			}

			//update upper diagonal part to 0
			if (j < i) {
				for (int k = 0; k < m; k++) {
					for (int p = 0; p < m; p++) {
						copy[j*m+k][i*m+p] = 0;
					}
				}
				continue;
			}

			//if j > i (the lower part of the matrix)

			//modify the rest of the ith column of this block of matrices
			//calculate the inverse of the diagonal chol factorization matrix
			vector<vector<double>> inverse_chol_blc(m, vector<double>(m));
			vector<vector<double>> transpose_chol_blc(m, vector<double>(m));
			inverse_chol_blc = lower_inverse(chol_block, m); 

			for (int k = 0; k < m; k++) {
				for (int p = 0; p < m; p++) {
					transpose_chol_blc[k][p] = inverse_chol_blc[p][k];
				}
			}

			//update the rest of blc mat along this column
			vector<vector<double>> blc(m, vector<double>(m));

			for (int k = 0; k < m; k++) {
				for (int p = 0; p < m; p++) {
					blc[k][p] = copy[j*m+k][i*m+p];//blc is A_{ji}
				}
			}
			blc = mat_mul(blc, transpose_chol_blc, m);
			for (int k = 0; k < m; k++) {
				for (int p = 0; p < m; p++) {
					copy[j*m+k][i*m+p] = blc[k][p];
				}
			}
		}
		//printf("after updating the ith block column\n");
		//printMatrix(copy, 4);//TODO

		//update the rest of the block matrices
		//k iterates through rest of blc mat in ith column
		for (int k = i+1; k <  num_block; k++) {
			vector<vector<double>> L_T(m, vector<double>(m));
			vector<vector<double>> tmp(m, vector<double>(m));

			for (int f = 0; f < m; f++) {
				for (int g = 0; g < m; g++) {
					tmp[f][g] = copy[k*m+f][i*m+g];
				}
			}
			//printMatrix(tmp,1);

			//p iterates through rest of blc mat in ith column
			// p can also be used to iterate through blc mat horizontally
			for (int p = i+1; p < num_block; p++) {
				vector<vector<double>> L(m, vector<double>(m));
				vector<vector<double>> tmp_2(m, vector<double>(m));
				for (int f = 0; f < m; f++) {
					for (int g = 0; g < m; g++) {
						L[f][g] = copy[p*m+f][i*m+g];
					}
				}
				//printMatrix(L,1);//TODO

				//transpose
				for (int f = 0; f < m; f++) {
					for (int g = 0; g < m; g++) {
						L_T[f][g] = L[g][f];
					}
				}

				//now update each blc mat left
				//can't directly store the matrix multiplication result
				//to tmp because inside a loop this tmp will not be updated
				//which results in the previous multiplication result being
				//multiplied again in a new iteration
				tmp_2 = mat_mul(tmp, L_T, m);
				for (int f = 0; f < m; f++) {
					for (int g = 0; g < m; g++) {
						copy[k*m+f][p*m+g] -= tmp_2[f][g];
					}
				}
				//printf("k: %d, p: %d.\n", k, p);
				//printf("updated the rest columns:\n");
				//printMatrix(copy, 4);//TODO

			}
		}
		//printf("%ith iteration in block_cho\n",i);//TODO
		//printMatrix(copy, 4);
		//printf("\n");//TODO

	}

	return copy;
}

/*
int main() {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> rand(0.0, 1.0);

	int size = 100;
	int m = 10;
	int num_block = floor(size/m);
	vector<vector<double>> A(size, vector<double>(size));
	vector<vector<double>> A_temp(size, vector<double>(size));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			A_temp[i][j] = rand(gen);

		}
	}
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			A[i][j] = 0.5 * (A_temp[i][j] + A_temp[j][i]);
		}
	}


	for (int i = 0; i < size; i++) {
		A[i][i] += 5;
	}

	//TODO: designate a matrix A
	//A = {{1,2,1,1},{2,13,8,8},{1,8,14,14},{1,8,14,15}};
	cout<<"original matrix: "<<endl;
	printMatrix(A, size);
	printf("\n");//TODO


	vector<vector<double>> alpha(size, vector<double>(size));
	alpha = block_chol(A, size, m, num_block);
	//alpha = cholesky(A, size);

	//cout<<"L*L^T: "<<endl;
	cout<<"cholesky_result L:"<<endl;//TODO

	vector<vector<double>> new_A (size, vector<double>(size));
	vector<vector<double>> alpha_T(size, vector<double>(size));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++){
			alpha_T[i][j] = alpha[j][i];
		}   
	}   
	printMatrix(alpha,size);
	printf("alpha_transpose: \n");//TODO
	printMatrix(alpha_T,size);
	new_A = mat_mul(alpha, alpha_T, size);
	printf("alpha times alpha_transpose\n");//TODO
	printMatrix(new_A, size);
	cout<<"norm: "<< norm(A, new_A, size)<<endl;
	//test1 for mat_mul
	vector<vector<double>> D = {{1,2,3,4},{4,3,2,1},{2,1,4,3},{3,4,1,2}};
	vector<vector<double>> B = {{1,2,1,2},{3,4,1,2},{2,1,2,3},{4,1,3,2}};
	vector<vector<double>> C(4, vector<double>(4));
	C = mat_mul(D,B, 4);
	printMatrix(C,4);
	//test2 for lower_inverse
	   vector<vector<double>> D = {{1.0/2,0,0,0},{-2.0/3,1.0/3,0,0},
	   {1.0/24,-1.0/12,1.0/4,0},{23.0/48,-11.0/24,-1.0/8,1.0/2}};
	   vector<vector<double>> B(4, vector<double>(4));
	   B = lower_inverse(D, 4);
	   printMatrix(B,4);
	return 0;
}
*/
