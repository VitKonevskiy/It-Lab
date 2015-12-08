#pragma once
#include "stdio.h"
#include <iostream>
#include <math.h>
using namespace std;
const double dt = 0.01;						//Еденица времени
const double eps = 0.01;					//
double a = 0;								//точка равновесия (по х)
const double D = 0.5;						//Связь
const int N = 30;							//Кол-во клеток
double I = 0.0;								//Импульс
#pragma warning(disable: 4996)
//Если у нас последовательность координат:
double f(double *x, double y, int i, double I)
{
	double c = 0;
	if (i == 0)
		c = D*(x[1] - x[0]);
	else if (i == N - 1)
		c = D*(x[N - 2] - x[N - 1]);
	else
		c = D*(x[i - 1] + x[i + 1] - 2 * x[i]);
	return (x[i] - (x[i] * x[i] * x[i] / 3) + I - y) + c;
}

/*
//Если у нас кольцо:
double f(double *x, double y, int i, double I)
{
double c = 0;
if (i == 0)
c = D*(x[N - 1] + x[1] - 2 * x[0]);
else if (i == N - 1)
c = D*(x[N - 2] + x[0] - 2 * x[N - 1]);
else
c = D*(x[i - 1] + x[i + 1] - 2 * x[i]);
return (x[i] - (x[i] * x[i] * x[i] / 3) + I - y) + c;//+и
}
*/

double g(double x, double y, int i)
{
	//Если меняем а:
	//return (eps*(x + (0.9 / N)*i));
	//Иначе:
	return (eps*(x + a));
}

int main()
{

	int q = 1;
	int p = 0;
	double T = 1000;
	double *x0 = new double[N];
	double *y0 = new double[N];
	double *x1 = new double[N];
	double *y1 = new double[N];
	double rkx[4][N];
	double rky[4][N];

	double changeA = -0.9 / T;				//Для Меняем а с течением времени

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < N; j++)
		{
			rkx[i][j] = 0;
			rky[i][j] = 0;
		}
	}
	for (int j = 0; j < N; j++)
	{
		x0[j] = 0.1 + 0.01*j;
		y0[j] = 0;
	}
	double time = 0;
	FILE *xfile = fopen("D:/It-lab/filexi2.txt", "w");
	while (time < T)
	{

		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < N; i++)
			{

				if (j == 3)
				{
					x1[i] = x0[i] + rkx[2][i];
					y1[i] = y0[i] + rky[2][i];
				}
				else if (j > 0)
				{
					x1[i] = x0[i] + rkx[2][i - 1] / 2;
					y1[i] = y0[i] + rky[2][i - 1] / 2;
				}
				else
				{
					x1[i] = x0[i];
					y1[i] = y0[i];
				}
			}
			for (int i = 0; i < N; i++) {
				if (i != 0) I = 0;		//Возбудили только первую клетку
				rkx[j][i] = dt*f(x1, y1[i], i, I);
				rky[j][i] = dt*g(x1[i], y1[i], i);
			}
		}
		for (int i = 0; i < N; i++)
		{
			//if (time == T / 4)			//Возвращаем значения возбужденных клеток
			//for (int k = N - 1; k>N / 4; k--)
			//	x0[k] = 0.01*k;
			x0[i] += (rkx[0][i] + 2 * rkx[1][i] + 2 * rkx[2][i] + rkx[3][i]) / 6;
			y0[i] += (rky[0][i] + 2 * rky[1][i] + 2 * rky[2][i] + rky[3][i]) / 6;
			//fprintf(xfile, " %g %g %g\n", time, x0[i], y0[i]);
			int itime = time * 100;
			if ((itime % 100) == 0)
			{
				fprintf(xfile, "%g ", x0[i]);
				if (i == 29)
				{
					fprintf(xfile, "\n");
					p++;
				}
			}
		}
		time += dt;
	}
	delete[] x0;
	delete[] y0;
	delete[] x1;
	delete[] y1;
	return 0;
};