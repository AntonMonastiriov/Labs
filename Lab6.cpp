// 4 варіант
// 
// 
#include <math.h>
#include <iostream>

using namespace std;
// tau 0.01
// a 0
// b 0.6
// c 0
// d 0.1
// h 0.05

double* ThreediagonalMatrix(double* A, double* B, double* C, double* D, int n)
{
	double* p = new double[n];
	double* q = new double[n];
	double* M = new double[n];
	p[0] = C[0] / B[0];
	q[0] = D[0] / B[0];
    for(int i = 1; i < n - 1; ++i)
    {
        p[i] = C[i] / (B[i] - A[i] * p[i - 1]);
        q[i] = (D[i] - A[i] * q[i - 1]) / (B[i] - A[i] * p[i - 1]);
    }
    q[n - 1] = (D[n - 1] - A[n - 1] * q[n - 2]) / (B[n - 1] - A[n - 1] * p[n - 2]);
    C[n - 1] = 0;
    M[n - 1] = q[n - 1];

    for(int i = n-2; i >= 0; --i)
    {
        M[i] = -p[i] * M[i + 1] + q[i];
    }
    return M;
}



double** ImplicitTwoLayer(double a, double b, double c, double d, double h, double tau)
{
    int N = (b - a) / h + 1;
    int M = (d - c) / tau + 1;
    double** U = new double*[N];
    for (int i = 0; i < N; i++) U[i] = new double[M];
    double x = a, t = c + tau;
    for(int i = 0; i < N; i++)
    {
        U[i][0] = sin(2*x);
        x += h;
    }
    for (int j = 1; j < M; j++)
    {
        U[0][j] = 2*t;
        U[N - 1][j] = 0.932;
        t += tau;
    }
    x = a + h; t = c + tau;
    double* aa = new double[N - 2];
    double* bb = new double[N - 2];
    double* cc = new double[N - 2];
    double* dd = new double[N - 2];
    for(int k = 0; k < N-2; k++)
    {
        aa[k] = tau / (h * h);
        bb[k] = -(1 + 2 * tau / (h*h));
        cc[k] = tau/ (h*h);
        dd[k] = -sin(2*x);
        x += h;
    }
    double* MM = ThreediagonalMatrix(aa, bb, cc, dd, N-2);

    for (int i = 1; i < N - 1; i++) U[i][1] = MM[i - 1];

    x = a + h; t = c + tau;
    for(int j = 2; j < M; j++)
    {
        x = a + h;
        for (int i = 1; i < N - 1; i++)
        {
            aa[i - 1] = tau / (h * h);
            bb[i - 1] = -(1 + 2 * tau / (h * h));
            cc[i - 1] = tau / (h * h);
            dd[i - 1] = -U[i][j - 1];
            x += h;
            MM = ThreediagonalMatrix(aa, bb, cc, dd, N-2);
            U[i][j] = MM[i - 1];
        }
        t += tau;
    }
    return U;
}

int main()
{
    double** U = ImplicitTwoLayer(0, 0.6, 0, 0.1, 0.05, 0.01);
    int N = 0.6 / 0.05 + 1;
    int M = 0.1 / 0.01 + 1;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M-1; j++)
        {
            printf("%.03f\t",U[j][i]);
            if (j == M - 1) printf("\n");
        }
    }

}
