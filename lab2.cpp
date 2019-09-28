/************************************************************************
*file: lab2.cpp
*author: Volodymyr Rekeda
*written: 05/12/2018
*last modified: 06/12/2018
*************************************************************************/

#include <iostream>
#include <cmath>
#include <windows.h>
#include <iomanip>

#define step 0.5

using namespace std;

//returns value of function for given x
double func(double x) {
	return x*x + 2 * sin(x) - 1;
}

//returns value of first derivative for given x
double func_first_derivative(double x) {
	return 2 * x + 2 * cos(x);
}

//returns value of second derivative for given x
double func_second_derivative(double x) {
	return 2 - 2 * sin(x);
}

//finds sub-interval which contains root at given interval (a; b)
double* find_roots_interval(double a, double b, int *n) {
	double *begs_of_intervals = nullptr;
	do {
		if (func(a)*func(a + step) < 0) {
			(*n)++;
			begs_of_intervals = (double*)realloc(begs_of_intervals, sizeof(double) * *n);
			begs_of_intervals[*n - 1] = a;
		}
		a += step;
	} while (a < b);
	return begs_of_intervals;
}

//calculates the root with given accurancy (eps) at given interval (a; a + step) using iteration method
double iteration_method(double a, int nmb, double eps, double *error, int* iter_nmb) {
	double M, m, lambda, q, xn, b;
	*iter_nmb = 0;
	b = a + step;
	double x0 = (a + b) / 2;

	if (fabs(func_first_derivative(a)) > fabs(func_first_derivative(b))){
		M = func_first_derivative(a);
		m = func_first_derivative(b);
	}
	else {
		M = func_first_derivative(b);
		m = func_first_derivative(a);
	}
	lambda = 1 / M;
	q = 1 - fabs(m / M);
	xn = x0;
	do {
		x0 = xn;
		xn = x0 - lambda * func(x0);
		(*iter_nmb)++;
	} while (fabs(xn - x0) > (1 - q) / q * eps);
	*error = (fabs(xn - x0))*q / (1 - q);
	return xn;
}

//calculates the root with given accurancy (eps) at given interval (a; a + step) using tangent method
double tangent_method(double a, int nmb, double eps, double *error, int *iter_nmb) {
	double m, x0, xn;
	double b = a + step;

	*iter_nmb = 0;
	m = fabs(func_first_derivative(a));

	for (double i = a; i < b; i += 0.001) {
		if (fabs(func_first_derivative(i) < m)) {
			m = fabs(func_first_derivative(i));
		}
	}
	if (func(a) * func_second_derivative(a) > 0)
		x0 = a;
	else
		x0 = b;
	xn = x0;
	do {
		xn = xn - func(xn) / func_first_derivative(xn);
		(*iter_nmb)++;
	} while (fabs(func(xn) / m) > eps);
	*error = fabs(func(xn)) / m;
	return xn;
}

void table_iter(double *begs_of_intervals, int nmb) {
	double err = 0, root;
	int iter_nmb = 0;

	cout << "\tIteration method:\n";
	cout << "----------------------------------------" << endl;
	cout << "|  eps |     root     |     error      |" << endl;
	cout << "----------------------------------------" << endl;

	for (int i = 0; i < nmb; i++) {
		for (double eps = 0.01; eps >= 1e-14; eps *= 1e-3) {
			root = iteration_method(begs_of_intervals[i], nmb, eps, &err, &iter_nmb);
			cout << "|" << setw(6) << eps 
				 << "|" << setprecision(10) << setw(14) << root
				 << "|" << setprecision(10) << setw(16) << err << "|" << endl;
			cout << "----------------------------------------" << endl;
		}
	}
	cout << endl;
}

void table_tangent(double *begs_of_intervals, int nmb) {
	double err = 0, root;	
	int iter_nmb = 0;

	cout << "\tTangent method:\n";
	cout << "----------------------------------------" << endl;
	cout << "|  eps |     root     |     error      |" << endl;
	cout << "----------------------------------------" << endl;

	for (int i = 0; i < nmb; i++) {
		for (double eps = 0.01; eps >= 1e-14; eps *= 1e-3) {
			root = tangent_method(begs_of_intervals[i], nmb, eps, &err, &iter_nmb);
			cout << "|" << setw(6) << eps
				 << "|" << setprecision(10) << setw(14) << root
				 << "|" << setprecision(10) << setw(16) << err << "|" << endl;
			cout << "----------------------------------------" << endl;
		}
	}
	cout << endl;
}

void table_compare(double *begs_of_intervals) {
	double err = 0;
	int iter_nmb1 = 0, iter_nmb2 = 0;

	cout << "\tNumber of iterations:\n";
	cout << "----------------------------------------" << endl;
	cout << "|  eps |Iteration method|Tangent method|" << endl;
	cout << "----------------------------------------" << endl;
	for (double eps = 0.01; eps >= 1e-14; eps *= 1e-3) {
		iteration_method(begs_of_intervals[0], 1, eps, &err, &iter_nmb1);
		tangent_method(begs_of_intervals[0], 1, eps, &err, &iter_nmb2);
		cout << "|" << setw(6) << eps
			 << "|" << setw(16) << iter_nmb1
			 << "|" << setw(14) << iter_nmb2 << "|" << endl;
		cout << "----------------------------------------" << endl;
	}
	cout << endl;
}

int main() {
	system("color F0");

	int nmb = 0;
	double *begs_of_intervals;
	double a = -10.0, b = 10.0;
	cout << "x*x + 2 * sin(x) - 1" << endl << endl;

	begs_of_intervals = find_roots_interval(a, b, &nmb);
	
	if (nmb > 0) {
		table_iter(begs_of_intervals, nmb);
		table_tangent(begs_of_intervals, nmb);
		table_compare(begs_of_intervals);
	}
	else
		cout << "No roots on interval: " << a << " " << b << endl;

	system("pause");
}