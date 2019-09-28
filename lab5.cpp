/************************************************************************
*file: lab4.cpp
*author: Volodymyr Rekeda
*written: 12/12/2018
*last modified: 16/12/2018
*var: 14 = 14 % 8 = 6 = (110)
*многочлени Чебишева
*схема з вибором головного елемента
*інтегрування коефіцієнтів нормальної системи за узагальненою формулою Сімпсона
*************************************************************************/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#define a 1//2
#define b 35//9
#define eps 10e-2
#define max_n 120

double *matrix_roots;

using namespace std;

double func(double x) {
	return log10(x*x)*sin(x / 2.0) * exp(pow(x, 1.0 / 7.0));
	//return 5.0 * log10(x)*sin(x)*cos(2.0 * x);
}

double* main_element(double **initial_matrix, const int n) {
	double max = 0;
	int p, q;
	double **M = new double*[n], **triangle_m = new double*[n];
	double *mi = new double[n], *roots = new double[n+1];
	int m = n + 1;

	for (int i = 0; i < n; i++) {
		M[i] = new double[m];
		triangle_m[i] = new double[m];
	}

	//init matrix M
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M[i][j] = initial_matrix[i][j];
			triangle_m[i][j] = 0;
		}
		M[i][n] = initial_matrix[i][max_n];
	}
	for (int i = 0; i <= n; i++)
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

	//cout << "Triangle matrix:\n";
	//print_matrix(triangle_m, n, m);

	//calculate roots and set them to their positions
	for (int i = n - 1; i >= 0; i--) {
		for (int j = 0; j < m - 1; j++)
			if ((triangle_m[i][j] != 0) && (roots[j] != 0))
				triangle_m[i][m - 1] -= triangle_m[i][j] * roots[j];

		for (int j = 0; j < m - 1; j++)
			if ((triangle_m[i][j] != 0) && (roots[j] == 0))
				roots[j] = triangle_m[i][m - 1] / triangle_m[i][j];
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

double phi(double x, int n) {
	if (n == 0)
		return 1.0;
	if (n == 1)
		return x;

	double t0 = 1, t1 = x;
	double buf = 0;
	for (int i = 1; i < n; i++)
	{
		buf = 2.0 * x * t1 - t0; ;
		t0 = t1;
		t1 = buf;
	}
	return t1;
}

double pm(double x, int m)
{
	double p = 0;
	for (int i = 0; i < m; i++)
		p += matrix_roots[i] * phi(x, i);

	return p;
}

double sigma(double x, int m, int)
{
	double result = func(x) - pm(x, m);

	return result * result;
}

double phi_2(double x, int i, int j) {
	if (i == j)
	{
		double res = phi(x, i);
		return res * res;
	}
	if (j == max_n)
	{
		return func(x) * phi(x, i);
	}
	return phi(x, i) * phi(x, j);
}

double simpson_method(double h, int n, int i, int j, double(*func)(double, int, int)) {
	
	double sigma_1 = 0.0, sigma_2 = 0.0;
	for (int k = 1; k < n; k += 2)
	{
		sigma_1 += func(a + k * h, i, j);
	}

	for (int k = 2; k < n - 1; k += 2)
	{
		sigma_2 += func(a + k * h, i, j);
	}

	return h / 3 * (func(a, i, j) + func(b, i, j) + 4 * sigma_1 + 2 * sigma_2);
}

double calc_integral(int i, int j, double(*func)(double, int, int)) {
	double eps_integr = 1e-3;
	int n = (b - a) / sqrt(sqrt(eps_integr));
	double h = (b - a) / (double)n;
	double integral_1n = 0, integral_2n = 0;

	n = (b - a) / h;
	if (n % 2 == 1)
		n++;
	h = (b - a) / (double)n;

	integral_1n = simpson_method(h, n, i, j, func);
	n *= 2;
	h /= 2;
	integral_2n = simpson_method(h, n, i, j, func); ;
	while (fabs((integral_1n - integral_2n) / integral_2n) > 15 * eps_integr)
	{
		integral_1n = integral_2n;
		n *= 2;
		h /= 2;
		integral_2n = simpson_method(h, n, i, j, func);
	}
	return integral_2n;
}

double **create_matrix(int n) {
	double **matrix = new double*[n];

	for (int i = 0; i < n; i++) 
		matrix[i] = new double[n + 1];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n + 1; j++)
			matrix[i][j] = calc_integral(i, j, phi_2);

	return matrix;
}

void print_matrix(double **matrix, int n) {

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n + 1; j++)
			cout<<matrix[i][j]<<" ";
		cout << endl;
	}
}

double get_delta(int m) {
	double h = (b - a) / 100.0;
	double sum = 0;
	for (double x = a; x <= b; x += h) {
		sum += pow((func(x) - pm(x, m)), 2);
	}
	return sqrt(sum / 101);
	//return sqrt(calc_integral(m, -1, sigma) / (b - a));
}

int main() {
	
	double **matrix;
	double delta;
	ofstream out_file("table.csv");
	cout << "calculating...\n";

	matrix = create_matrix(max_n);

	for (int i = 1; i < max_n; i++) {
		matrix_roots = main_element(matrix, i);
		delta = get_delta(i + 1);

		cout << "delta = " << delta << " m = " << i << endl;
		if (delta < eps) {
			double h = (b - a) / 100.0;
			double x = a;
			cout << "Least Squares Deviation = " << delta << " for m = " << i << endl;
			cout << "  x\tpm\n";
			cout.setf(ios_base::fixed);
			out_file.setf(ios_base::fixed);

			for (int j = 0; j <= 100; j++) {
				x = a + j * h;
				cout << setprecision(2) << x << " " << setprecision(8) << pm(x, i + 1) << endl;
				out_file << setprecision(3) << x << ";" << setprecision(8) << pm(x, i + 1) << endl;
			}

			delete[] matrix_roots;
			break;
		}
		delete[] matrix_roots;
	}

	for (int i = 0; i < max_n; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
	out_file.close();
	system("pause");
	return 0;
}
