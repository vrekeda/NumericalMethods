/************************************************************************
*file: lab6.cpp
*author: Volodymyr Rekeda
*written: 17/12/2018
*last modified: 21/12/2018
*var: 14
*************************************************************************/

#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

double acc_solution(double x) {
	return x * x*x + 3 * x + 1;
}

double f(double x, double y, double z) {
	return 2 * x*z / (x*x + 1);
}

double g(double x, double y, double z) {
	return z;
}

double Runge_Kutta_4(double h, double a, double b) {
	double y = 1.0, z = 3.0;
	double k1, k2, k3, k4;
	double q1, q2, q3, q4;
	double x;
	for (x = a; x < b;) {
		k1 = h * f(x, y, z);
		q1 = h * g(x, y, z);

		k2 = h * f(x + h / 2.0, y + q1 / 2.0, z + k1 / 2.0);
		q2 = h * g(x + h / 2.0, y + q1 / 2.0, z + k1 / 2.0);

		k3 = h * f(x + h / 2.0, y + q2 / 2.0, z + k2 / 2.0);
		q3 = h * g(x + h / 2.0, y + q2 / 2.0, z + k2 / 2.0);

		k4 = h * f(x + h, y + q3, z + k3);
		q4 = h * g(x + h, y + q3, z + k3);
		x += h;
		z = z + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
		y = y + (q1 + 2.0*q2 + 2.0*q3 + q4) / 6.0;
	}
	return y;
}

double task_1(double *h, double eps, double a, double b) {
	double accurate_solution = acc_solution(b);
	double delta;
	*h = (b - a) / 10.0;
	delta = fabs(accurate_solution - Runge_Kutta_4(*h, a, b));
	while (delta > eps) {
		*h /= 2;
		delta = fabs(accurate_solution - Runge_Kutta_4(*h, a, b));
	}
	//cout << delta << " " << Runge_Kutta_4(*h, a, b)<<" "<<*h;
	return delta;
}

double task_2(double eps, double a, double b) {
	double h = sqrt(sqrt(eps));
	double y_1, y_2;
	y_1 = Runge_Kutta_4(h, a, b);
	h /= 2;
	y_2 = Runge_Kutta_4(h, a, b);
	while (fabs(y_1 - y_2) > 15 * eps) {
		y_1 = y_2;
		h /= 2;
		y_2 = Runge_Kutta_4(h, a, b);
	}
	return h;
}

void task_3(double a, double b) {
	double n = 10;
	double h = (b - a) / n;
	cout.setf(std::ios_base::fixed);
	cout << "  x\t   y\n";
	for (double x = a; x <= b; x += h) {
		cout << setprecision(2) << x << "\t" << setprecision(6) << Runge_Kutta_4(h, a, x) << endl;
	}
}

int main() {
	double a = 0, b = 1;
	double h, h_2, eps_1 = 1e-2, eps_2;

	eps_2 = task_1(&h, eps_1, a, b);
	cout << "Task 1:\n";
	cout << "val = " << Runge_Kutta_4(h, a, b) << endl;
	cout << "delta = " << eps_2 << endl << endl;

	h_2 = task_2(eps_2, a, b);
	cout << "Task 2:\n";
	cout << "h in 1st task - " << h << "\nh in 2nd task - " << h_2 << endl;

	cout << "\nTask 3:\n";
	task_3(a, b);

	system("pause");
	return 0;
}