// Варіант 10

// y' = y*(ln(y/x)+1)/x |x0 = 1|y0 = 1.1|xn = 7|y = x*exp(Cx)
#include <iostream>
#include <math.h>
using namespace std;
double* Runga_Kutta();
double Af(double x);
double f(double x, double y);

#define h 0.1
#define a 1.0
#define b 7.0
int main()
{
	Runga_Kutta();
}


double f(double x, double y)
{
	return y * (log(y / x) + 1) / x;
}

double Af(double x)
{
	return x * exp(x * log(1.1));
}

double* Runga_Kutta()
{
	int n = (b - a) / h + 1;
	double* x = new double[n];
	double* y = new double[n];
	double* k1 = new double[n];
	double* k2 = new double[n];
	double* k3 = new double[n];
	double* k4 = new double[n];
	double* d = new double[n];
	x[0] = 1;
	y[0] = 1.1;
	printf("Runge-Kutta 4th order\nx\tF(x)\t\tAprox\n");
	for (int i = 1; i < n; i++)
	{
		k1[i] = h * f(x[i-1], y[i-1]);
		k2[i] = h * f(x[i-1] + h / 2, y[i-1] + k1[i] / 2);
		k3[i] = h * f(x[i-1] + h / 2, y[i-1] + k2[i] / 2);
		k4[i] = h * f(x[i-1] + h, y[i-1] + k3[i]);
		y[i] = y[i - 1] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
		x[i] = x[i - 1] + h;
		printf("%.2f\t%.6f\t%.6f\n", x[i], Af(x[i]), y[i]);
	}
	return y;
}