/*  Methods Used For Solving Numerical Integral  */

#ifndef _QUADRULES_H
#define _QUADRULES_H

#define _USE_MATH_DEFINES 
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>

double cheb(int n, double x);
double evalLeg(int n, double x);
void GauLeg(int n, double* x);
void swap(double& a, double& b);
void gauss_jordan(double* a, double* b, int n, int m);
void gauss_jordan(double* a, int n);
void printMatrix(double* a, int m, int n);
void multMatrix(double* a, double* b, double* erg, int n_row, int n_sum, int n_col);


template <class T>
double trapz(T& func, const double a, const double b, const int n )
{
    double h = (b - a) / n;
    double s = 0;
    double sum = 0;
    for (int i = 1; i < n; i++)  s += func(a + i * h);
    sum += h * (func(a) + 2 * s + func(b)) / 2;
    return sum;
}

template <class T>
double trapz(T& func, const double xi, const double a, const double b, const int n)
{
    double h = (b - a) / n;
    int n1 = (int)((xi - a) / h + 0.5);
    int n2 = n1 + 1;
    double s1 = 0;
    double s2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    double sum = 0;
    if (n1 != 0 && n1 != n)
    {
        for (int i = 0; i < n1; i++)
        {
            s1 += func(a + i * h);
        }
        sum1 = h * s1 - 0.5 * h * func(a) + 0.5 * h * func(xi - h);
        for (int i = n2; i <= n; i++)
        {
            s2 += func(a + i * h);
        }
        sum2 = h * s2 - 0.5 * h * func(b) + 0.5 * h * func(xi + h);
        sum = sum1 + sum2;
    }
    else if (n1 == 0)
    {
        for (int i = 1; i <= n; i++)
        {
            s1 += func(a + i * h);
        }
        sum = h * s1 - 0.5 * h * func(b) + 0.5 * h * func(a + h);
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            s1 += func(a + i * h);
        }
        sum = h * s1 - 0.5 * h * func(a) + 0.5 * h * func(b - h);
    }

    return sum;

}


template <class T>
double gausslegend(T& func, const double a, const double b, const int n)
{
    double sum = 0;
    double* x = new double[2*n];
    double c1 = (b - a) / 2;
    double c2 = (b + a) / 2;
    GauLeg(n, x);
    if (a != b)
    {
        for (int i = 0; i < n; i++)
        {
            sum += c1 * x[2 * i + 1] * func(c2 + c1 * x[2 * i]);
        }
    }
    else
    {
        sum = 0;
    }
    
    delete [] x;
    return sum;
}























#endif
