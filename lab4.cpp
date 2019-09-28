/************************************************************************
*file: lab4.cpp
*author: Volodymyr Rekeda
*written: 11/12/2018
*last modified: 12/12/2018
*var: 14
*************************************************************************/

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

#define b 27
#define a 1

double func(double x) {
	return sin(x)*sin(x);
}

double func_primitive(double x) {
	return -sin(2 * x) / 4 + x / 2;
}

//calculates value of integral using trapezium method
double trapezium_method(int n) {
	double result = func(a) / 2;
	double h = (b - a) / (double)n;
	for (int i = 1; i < n; i++) {
		result += func(a + h * i);
	}
	result += func(b) / 2;
	return result * h;
}

double table1(double integral0) {
	double delta;
	double eps = 0.0001;
	double h;
	int n;
	double integral_n;

	h = sqrt(12 * eps / (b - a) / 2);//calculate h, /2 because 2nd derivative of function is 2(cos^2(x) - sin^2(x)) 
	n = (b - a) / h;
	integral_n = trapezium_method(n);
	delta = fabs(integral0 - integral_n);

	cout << "Task 1" << endl;
	cout << " eps  |      h      |    I     |     delta" << endl;
	cout<< eps 
		<< "|" << (b - a) / (double)n 
		<< "|" << integral_n 
		<< "|" << delta << endl;
	return delta;
}

void table2(double delta, double integral) {
	int n;
	double integral_2n, integral_1n;

	n = (b - a) / sqrt(delta);
	integral_1n = trapezium_method(n);
	n *= 2;
	integral_2n = trapezium_method(n);
	while (fabs(integral_2n - integral_1n) / 3 > delta) {
		integral_1n = integral_2n;
		n *= 2;
		integral_2n = trapezium_method(n);
	}

	cout << "\nTask 2" << endl;
	cout << "     delta     |       h       |   eps" << endl;
	cout << delta
		<< "|" << (b - a) / (double)n
		<< "|" << fabs(integral - integral_2n) << endl;
}

int main() {
	double integral = func_primitive(b) - func_primitive(a); //calculate value of integral
	double delta;
	cout << "Integral from 1 to 27 sin^2(x) dx\nF = -sin(2x)/4 + x/2\n";
	cout << "I = " << setprecision(10) << integral << endl;
	delta = table1(integral);
	table2(delta, integral);
	system("pause");
	return 0;
}