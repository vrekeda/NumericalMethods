/************************************************************************
*file: lab1.cpp
*author: Volodymyr Rekeda
*written: 25/09/2018
*last modified: 26/09/2018
*************************************************************************/

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

#define a -31.3
#define b -4.9

int task1_cos(long double x, long double eps, long double  *sum, long double *u){
	double n = 0, i = 0;
	x = fabs(x);
	while (x > 2 * M_PI) {
		x -= 2 * M_PI;
	}
	*sum = 0;
	*u = 1;
	while (eps < fabs(*u)) {
		*sum += *u;
		*u = (-1) * *u * x * x / (i + 1) / (i + 2);
		i += 2;
		n++;
	}
	return n;
}

void print_task1() {
	long double rest = 0, result = 0;
	int n;
	long double x = (b + a) / 2;
	std::cout << "   eps  |  n |    delta    |     rest" << std::endl;
	for (long double eps = 1e-2; eps >= 1e-14; eps *= 1e-3) {
		n = task1_cos(x, eps, &result, &rest);
		std::cout.width(7);
		std::cout << eps << " | ";
		std::cout.width(2);
		std::cout << n << " | ";
		std::cout.width(11);
		std::cout << fabs(result - cos(x)) << " | " << rest << std::endl;
	}
}

long double task2_cos(int n, long double x, long double *u) {
	long double result = 0;
	long double k = 0;
	*u = 1;
	x = fabs(x);
	while (x > 2 * M_PI) { 
		x -= 2 * M_PI;
	}
	for (int i = 0; i < n; i++) {
		result += *u;
		*u = (-1) * *u * x * x / (k + 1) / (k + 2);
		k += 2;
	}
	return result;
}

void print_task2() {
	long double rest = 0, result = 0, x = a;
	int n;
	long double h = (b - a) / 10;
	n = task1_cos((b + a) / 2, 1e-8, &result, &rest);
	std::cout << "n=" << n << std::endl;
	std::cout << "     x  |    delta    |     rest" << std::endl;
	for (int i = 0; i <= 10; i++) {
		result = task2_cos(n, x, &rest);
		std::cout.width(7);
		std::cout << x << " | ";
		std::cout.width(11); 
		std::cout << fabs(result - cos(x)) << " | " << rest << std::endl;
		x += h;
	}
}

int main() {
	system("color F0");
	std::cout << "Task 1\n";
	print_task1();
	std::cout << "\nTask 2\n";
	print_task2();
	system("pause");
	return 0;
}

