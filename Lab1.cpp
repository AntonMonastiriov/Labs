// Перша лаба - варіант 6, до 17 f(x)=ln(x^4 - 2x^2+3) [a,b] = [-1.5, 1.5]; h= 0.12, linear
#include <cstdlib>
#include <iostream>
#include <math.h>
using namespace std;
int main()
{
    double x[9] = { -1.6,-1.1,-0.6,-0.2,0.2,0.7,0.93,1.2,1.51 };
    double y[9] = {1.48921,0.71496, 0.87946, 1.07213, 1.07213, 0.81541, 0.70223, 0.78554, 1.29161};
    double* X = new double[26];
    for (int i = 0; i < 26; i++)
    {
        X[i] = -1.5 + 0.12 * i;
    }
    double* a = new double[8];
    double* b = new double[8];
    double* trueY = new double[26];
    double* approxY = new double[26];
    double* Ly = new double[26];
    // Значення функції в точках 
    for (int i = 0; i <= 25; i++)
    {
        trueY[i] = log(X[i] * X[i] * X[i] * X[i] - 2 * X[i] * X[i] + 3);
    }
    // Кускова лінійна інтерполяція:
    for (int i = 0; i < 8; i++)
    {
        a[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        b[i] = y[i] - a[i] * x[i];
    }
    int j = 0;
    for (int i = 0; i < 8; i++)
    {
        for (; j <= 25; j++)
        {
            if (X[j] >= x[i] && X[j] <= x[i + 1])
            {
                approxY[j] = a[i] * X[j] + b[i];
            }
            else break;
        }
    }
    // Інтерполяція поліномом Лагранжа
    double multiplyer = 1;
    for (int k = 0; k <= 25; k++)
    {
        Ly[k] = 0;
        for (int i = 0; i < 9; i++) 
        {
            double multiplyer = 1;
            for (int j = 0; j < 9; j ++)
            {
                if (i != j)
                {
                    multiplyer *= (X[k] - x[j]);
                    multiplyer /= (x[i] - x[j]);
                }
            }
            Ly[k] += trueY[i] * multiplyer;
        }
    }

    printf("    x ------ y ------- apr.Y ---- Lagrang\n");
    for (int i = 0; i <= 25; i++)
    {
        printf("%.02f --- %.05f --- %.05f --- %.05f\n", X[i], trueY[i], approxY[i], Ly[i]);
    }

    // Диференціювання
    // 4(x*x*x)/(x^4+2x^2+3)
    double *df = new double[26];
    double* Apdf = new double[26];
    double* d2f = new double[26];
    double* Apd2f = new double[26];
    // Значення похідної, обчислене за формулами
    for (int i = 0; i <= 25; i++)
    {
        df[i] = 4 * (X[i] * X[i] * X[i]) / (X[i] * X[i] * X[i] * X[i] + 2 * X[i] * X[i] + 3);
        d2f[i] = (-4 * (3 - 7 * X[i] * X[i] - X[i] * X[i] * X[i] * X[i] + X[i] * X[i] * X[i] * X[i] 
            * X[i] * X[i])) / ((3 - 2 * X[i] * X[i] * +X[i] * X[i] * X[i] 
                * X[i]) * (3 - 2 * X[i] * X[i] * +X[i] * X[i] * X[i] * X[i]));
    }
    // Значення похідної, обчислене чисельно за допомогою кускової інтерполяції
    for (int i = 0; i <= 25; i++)
    {
        if (i == 0) Apdf[0] = (approxY[1] - approxY[0]) / 0.12;
        else if (i == 25) Apdf[25] = (approxY[25] - approxY[24]) / 0.12;
        else Apdf[i] = (approxY[i + 1] - approxY[i - 1]) / (0.24);
    }
    Apd2f[0] = 0;
    Apd2f[25] = 0;
    for (int i = 1; i < 25; i++)
    {
        Apd2f[i] = (approxY[i + 1] + approxY[i - 1] - 2*approxY[i]) / 0.0144;
    }

    printf("   x ----- y' -------- Approx y' \n");
    for (int i = 0; i <= 25; i++)
    {
        printf("%.02f --- %.05f --- %.05f\n", X[i], df[i], Apdf[i]);
    }

    printf("   x ----- y'' -------- Approx y''\n");
    for (int i = 0; i <= 25; i++)
    {
        printf("%.02f --- % .05f --- % .05f\n", X[i], d2f[i], Apd2f[i]);
    }



    delete[] a;
    delete[] b;
    delete[] trueY;
    delete[] approxY;
    delete[] Ly;
    delete[] X;
    delete[] df;
    delete[] d2f;
    delete[] Apdf;
    delete[] Apd2f;
}