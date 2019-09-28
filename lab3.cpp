/************************************************************************
*file: lab3.cpp
*author: Volodymyr Rekeda
*written: 07/12/2018
*last modified: 08/12/2018
*var: 14
*************************************************************************/

#include <iostream>
#include <Windows.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

void print_matrix(double *matrix, int n, int m) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			cout << setw(4) << matrix[i * m + j] << " ";
		cout << endl;
	}
	cout << endl;
}

void print_matrix(double **matrix, int n, int m) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			cout << setprecision(5) << setw(7) << matrix[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}


double* main_element(double *initial_matrix, const int n, const int m) {
	double max = 0;
	int p, q;
	double **M = new double*[n], **triangle_m = new double*[n];
	double *mi = new double[n], *roots = new double[n];


	for (int i = 0; i < n; i++) {
		M[i] = new double[m];
		triangle_m[i] = new double[m];
	}

	//init matrix M
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++) {
			M[i][j] = initial_matrix[i*m + j];
			triangle_m[i][j] = 0;
		}

	for (int i = 0; i < n; i++)
		roots[i] = 0;

	for (int k = 0; k < n; k++) {
		max = 0;

		//find main(max) element of matrix
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m - 1; j++)
				if (fabs(M[i][j]) > fabs(max)) {
					max = M[i][j];
					p = i;
					q = j;
				}

		//calculate multipliers
		for (int i = 0; i < n; i++)
			mi[i] = -M[i][q] / max;

		//remember main row
		for (int j = 0; j < m; j++)
			triangle_m[k][j] = M[p][j];


		for (int i = 0; i < n; i++)
			if (i != p)
				for (int j = 0; j < m; j++)
					M[i][j] += M[p][j] * mi[i];

		//set 0 to row and column which contain main element
		for (int i = 0; i < n; i++)
			M[i][q] = 0;
		for (int j = 0; j < m; j++)
			M[p][j] = 0;

	}

	cout << "Triangle matrix:\n";
	print_matrix(triangle_m, n, m);

	//calculate roots and set them to their positions
	int r_nmb;
	for (int i = n - 1; i >= 0; i--) {
		for (int j = 0; j < m - 1; j++) {
			if (roots[j] != 0) 
				triangle_m[i][m - 1] -= triangle_m[i][j] * roots[j];
			else if (triangle_m[i][j] != 0) 
				r_nmb = j;
		}
		roots[r_nmb] = triangle_m[i][m - 1] / triangle_m[i][r_nmb];
	}

	delete[] mi;
	for (int i = 0; i < n; i++) {
		delete[] M[i];
		delete[]triangle_m[i];
	}
	delete[] M;
	delete[] triangle_m;
	return roots;
}

double* simple_iteration(double *initial_matrix, const int n, const int m, int *k, double eps) {
	double *b = new double[n];
	double *xk = new double[m - 1];
	double *x = new double[m - 1];
	double **a = new double*[n];
	double q = 0, sum = 0, m_norm = 0;
	for (int i = 0; i < n; i++)
		a[i] = new double[m];

	//calculate vector b and x0
	for (int i = 0; i < n; i++)
		x[i] = b[i] =  initial_matrix[i*m + m - 1] / initial_matrix[i*m + i] ;

	//fill matrix of coefficients
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m-1; j++)
			if (i != j)
				a[i][j] = -initial_matrix[i*m + j] / initial_matrix[i*m + i];
			else
				a[i][j] = 0;

	//choose q
	for (int i = 0; i < n; i++)	{
		sum = 0;
		for (int j = 0; j < m-1; j++)
			sum += fabs(a[i][j]);
		if (sum > q)
			m_norm = sum;
	}
	q = m_norm;

	*k = 1;

	while (1) {
		
		m_norm = 0;

		//calculate x[k]
		for (int i = 0; i < n; i++) {
			sum = 0;
			for (int j = 0; j < m-1; j++){
				if (i != j)
					sum += a[i][j] * x[j];
			}
			xk[i] = b[i] + sum;
		}

		//calculate m-norm of vector |x[k]-x[k-1]|
		for (int i = 0; i < n; i++) {
			sum = 0;
			for (int j = 0; j < m-1; j++){
				sum += fabs(xk[j] - x[j]);
			}
			if (sum > m_norm)
				m_norm = sum;
			
		}

		//check accurancy
		if (m_norm <= eps * (1 - q) / q)
			break;

		for (int i = 0; i < n; i++) {
			x[i] = xk[i];
		}
		(*k)++;
	}

	delete[] b;
	delete[] x;
	for (int i = 0; i < n; i++)
		delete[] a[i];
	delete[] a;
	return xk;
}

int main() {
	const int N = 4;
	const int M = 5;
	int k = 0;
	double eps;
	double initial_matrix[N][M] = { {3,  19, 11, 8,  149},
									{9,  31, 3,  18, 257},
									{11, 7,  32, 13, 143},
									{12, 19, 12, 5,  144} };

	double modified_matrix[N][M] = { {-9, 0, -1, 3, 5 },
									 {9, 31, 3, 18, 257 },
									 {11, 7, 32, 13, 143 },
									 {12, 19, 12, 45, 344} };

	double *roots1, *roots2;

	cout << "Scheme with the choice of the main element:\n";
	cout << "Initial matrix:\n";
	print_matrix(initial_matrix[0], N, M);
	roots1 = main_element(initial_matrix[0], N, M);
	cout << "Roots:\n";
	for (int i = 0; i < N; i++) {
		cout << roots1[i] << " ";
	}

	eps = 1e-03;
	cout << "\n\nMethod of simple iterations:\n";
	cout << "Modified matrix:\n";
	print_matrix(modified_matrix[0], N, M);
	roots2 = simple_iteration(modified_matrix[0], N, M, &k, eps);
	cout << "Roots with accurancy eps = " << eps << ":\n";
	for (int i = 0; i < N; i++) {
		cout <<setprecision(7)<< roots2[i] << " ";
	}
	cout << endl;
	cout << "Number of iterations: " << k << endl << endl;

	delete[] roots1;
	delete[] roots2;

	system("pause");
	return 0;
}