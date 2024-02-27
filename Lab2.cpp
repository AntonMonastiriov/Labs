// Друга лаба - варіант 11, до 24
// Integrate[Sin(ax)/(a+x^2),{x,0,∞}], a = 1;2;...;10, h = 0.1, acc = 10^(-4), Simpson
#include <iostream>
#include <cstdlib>
#include <math.h>
using namespace std;
double Simpson(double h, double i, int j, double x);

int main()
{
	double x;
	double b = 20000;
	double acc = 0.0001;
	double* Fx = new double[10];
	double h = 0.01;
	int a = 10;
		for (int j = 1; j <= a; j++)
		{
			do 
			{
				Fx[j - 1] = 0;
				for (double i = 0.0; i <= b/h - 1; i++)
				{
					x = i * h;
					Fx[j - 1] += Simpson(h,i,j,x);
				}
				Fx[j - 1] *= h / 6;
				b += 1;
			} while (abs(h/6* Simpson(h, 1, j, b) - h/6* Simpson(h, 1, j, b+1))>acc);
		}
	printf(" a ----- F(x)\n");
	for (int i = 0; i < a; i++)
	{
		printf(" %d\t%.04f\n",i+1, Fx[i]);
	}
	delete[] Fx;
}
double Simpson(double h, double i, int j, double x)
{
	return sin(j*x)/(j+x*x)+4*(sin(j*(x+x+i*h)/2)/(j+((x+x+i*h)*(x+x+i*h))/4)+sin(j*(x+i*h))/(j+(x+i*h)*(x+i*h)));
}