// 4 варіант 
#include <cstdlib>
#include <math.h>
#include <iostream>
using namespace std;

double phi_y(double x);
double phi_x(double y);
void MPI(double& x0, double& y0);
void Seidel(double& x0, double& y0);
double* Jacobi(double A[][3], double* b);
void SeidelS(double A[][3], double* b);
int main()
{
	// Трансцендента система, метод простих іперацій
	// y = 1.5 - cos(x)
	// x = (1+sin(y-0.5))/2
	double x0 = 1, y0 = 1;
	int i = 0;
	MPI(x0, y0);
	double f1 = 1.5 - cos(x0) - y0;
	double f2 = (1 + sin(y0 - 0.5)) / 2 - x0;
	printf("x = %.04f\ty = %.04f\n", x0, y0);
	printf("Residual:\nf1 = %.08f\tf2 = %.08f\n", f1, f2);

	// Метод Зейделя
	x0 = 1, y0 = 1;
	Seidel(x0, y0);
	f1 = 1.5 - cos(x0) - y0;
	f2 = (1 + sin(y0 - 0.5)) / 2 - x0;
	printf("x = %.04f\ty = %.04f\n", x0, y0);
	printf("Residual:\nf1 = %.08f\tf2 = %.08f\n\n", f1, f2);

	// СЛАР
	double freeArr[3] = { 9.8,6.7,5.8 };
	double kArr[3][3]
	{
		{9.1,5.6,7.8},
		{4.1,5.7,1.2},
		{3.8,5.1,2.8}
	};

	double* x = Jacobi(kArr, freeArr);
	for (int i = 0; i < 3; i++) 
	{
		printf("x%i = %.04f\n", i+1,x[i]);
	}
	f1 = freeArr[0] - (kArr[0][0] * x[0] + kArr[0][1] * x[1] + kArr[0][2] * x[2]);
	f2 = freeArr[1] - (kArr[1][0] * x[0] + kArr[1][1] * x[1] + kArr[1][2] * x[2]);
	double f3 = freeArr[2] - (kArr[2][0] * x[0] + kArr[2][1] * x[1] + kArr[2][2] * x[2]);
	printf("Residual:\nf1 = %.08f\tf2 = %.08f\tf3 = %.08f\n\n", f1, f2, f3);

	SeidelS(kArr, freeArr);	
}

double phi_y(double x) 
{
	return 1.5 - cos(x);
}

double phi_x(double y)
{
	return (1+sin(y-0.5))/2;
}

void MPI(double& x0, double& y0)
{
	cout << "MPI\n";
	int i = 0;
	double x;
	double y;
	printf("k\txk\tyk\t|dk|\n");
	do
	{
		x = phi_x(y0);
		y = phi_y(x0);
		printf("%i\t%.04f\t%.04f\t%.04f\n", i + 1, x, y, max(abs(x - x0), abs(y - y0)));
		if (max(abs(x - x0), abs(y - y0)) < 0.0001) break;
		x0 = x;
		y0 = y;
		i++;
	} while (i < 1000);
}

void Seidel(double& x0, double& y0)
{
	cout << "Seidel\n";
	int i = 0;
	double x;
	double y;
	printf("k\txk\tyk\t|dk|\n");
	do
	{
		x = phi_x(y0);
		y = phi_y(x);
		printf("%i\t%.04f\t%.04f\t%.04f\n", i + 1, x, y, max(abs(x - x0), abs(y - y0)));
		if (max(abs(x - x0), abs(y - y0)) < 0.0001) break;
		x0 = x;
		y0 = y;
		i++;
	} while (i < 1000);
}

double* Jacobi(double A[][3], double* b)
{
	int k = 0;
	double x[3] = { 0.8, 0.5, -0.1 };
	double* xOld = new double[3];
	double* dx = new double[3];
	double norm;
	printf("Jacobi:\nk\tx1\tx2\tx3\tdk\n");
	do
	{
		k++;
		for (int i = 0; i < 3; i++)
		{
			xOld[i] = x[i];
			x[i] = b[i];
			for (int j = 0; j < 3; j++)
			{
				if (j != i)
				{
					x[i] -= A[i][j] * x[j];
				}
			}
			x[i] = x[i] / A[i][i];
			dx[i] = x[i] - xOld[i];
		}
		norm = max(abs(dx[0]),max(abs(dx[1]),abs(dx[2])));
		printf("%i\t%.04f\t%.04f\t%.04f\t%.04f\n", k, x[0], x[1], x[2], norm);
	} while (norm > 0.0001 && k < 1000);
	return x;
}

void SeidelS(double A[][3], double* b)
{
	int k = 0;
	double x[3] = { 0.7, 0.2, -1 };
	double* xOld = new double[3];
	double* dx = new double[3];
	double norm;
	printf("Seidel Slar:\nk\tx1\tx2\tx3\tdk\n");
	do
	{
		for (int i = 0; i < 3; i++)
		{
			xOld[i] = x[i];
		}
		k++;
		for (int i = 0; i < 3; i++)
		{
			x[i] = b[i];
			for (int j = 0; j < i; j++)
			{
				x[i] -= A[i][j] * x[j];
			}
			for (int j = i + 1; j < 3; j++)
			{
				x[i] -= A[i][j] * xOld[j];
			}
			x[i] = x[i] / A[i][i];
			dx[i] = x[i] - xOld[i];
		}
		norm = max(abs(dx[0]), max(abs(dx[1]), abs(dx[2])));
		printf("%i\t%.04f\t%.04f\t%.04f\t%.04f\n", k, x[0], x[1], x[2], norm);
	} while (norm > 0.0001 && k < 1000);
	for (int i = 0; i < 3; i++) printf("x%i = %.04f\t", i + 1, x[i]);
	double f1 = b[0] - (A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2]);
	double f2 = b[1] - (A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2]);
	double f3 = b[2] - (A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2]);
	printf("\nResidual:\nf1 = %.08f\tf2 = %.08f\tf3 = %.08f\n\n", f1, f2, f3);
}