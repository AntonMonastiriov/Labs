#include <iostream>
#include <cstdlib>
#include <math.h>
using namespace std;
double findXbisection(double& a, double& b);
double f(double x);
double df(double x);
double d2f(double x);
double findXNewton(double& x);
double findXNewtonb2(double x);
// Варіант 7
// e^x - e^(-3x) - 4 = 0
// f = e^x + e^(-3x) - 4
// f''= e^x + 9e^(-3x)
int main()
{
	cout << "Bisection method" << endl;
	double a = 0;
	double b = 2;
	double x = 1;
	int k = 1;

	printf("k ------ a ------- b ------ |b-a| ------ f(x)\n");
	do
	{
		x = findXbisection(a, b);
		printf("%d --- %.04f --- %.04f --- %.04f --- %.04f\n", k, a, b, abs(b-a), f(x));
		k++;
	} while (abs(f(x))>0.0001);
	cout << x << "\t" << f(x) << "\n";

	a = 0;
	b = 2;
	if (f(a) * d2f(a) > 0) x = a;
	else x = b;
	int i = 0;

	cout << endl << "Newton method" << endl;

	cout << "x(k) - x(k-1)\n";
	do
	{
		x = findXNewton(x);
		cout << abs(findXNewtonb2(x) - findXNewtonb2(x + (f(x)) / df(x))) << endl;
		i++;
	} while (abs(findXNewtonb2(x) - findXNewtonb2(x + (f(x)) / (df(x)))) > 0.0001);
	
	cout << x << "\t" << f(x) << "\n";
}

double f(double x)
{
	return exp(x) + exp(-3 * x) - 4;
}

double df(double x)
{
	return exp(x) - 3*exp(-3 * x);
}

double d2f(double x)
{
	return exp(x) + 9*exp(-3 * x);
}

double findXbisection(double &a, double &b)
{
	double x = (a+b)/2;
	if (f(x) == 0) return x;
	else if (f(a) * f(x) < 0) b = x;
	else a = x;
}

double findXNewton(double& x)
{	
	x = x - (f(x)) / (df(x));
	return x;
}
double findXNewtonb2(double x)
{
	return x - (f(x)) / (df(x));
}