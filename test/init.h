/*  initial parameters and functions  */

#ifndef _INIT_H
#define _INIT_H

#define _USE_MATH_DEFINES 
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <ctime>
#include <iomanip>
#include<omp.h>

#include "quadrules.h"

using namespace std;
using namespace Eigen;

/*  Initialize The Parameters  */
const int n = 50;
const double a = 0;
const double b = M_PI;
const double aa = 1e-1;
const double A = 1e10;
const double a1 = 2.086;
const double nf = 1;
const double c = log(10) * 3 * a1 / (2 * M_PI * M_PI);
const double c1 = nf * a1 / (3 * M_PI);
const double h1 = (A - aa) / n;
const double t1 = -1;
const double t2 = 10;
const double h2 = (t2 - t1) / n;

/*  All The Functions Or Variables That Will Be Used  */
double f(double, double, double(*fx)(VectorXd, double), VectorXd, double);
double f2(double);
double f3(double);
double f4(double);
double f5(double);
double pa(int, double);
double savx;
double Sy;
int I;
double (*intfunc)(double, double, double(*fx)(VectorXd, double), VectorXd, double);
double (*intfunc1)(double, double, double(*fx)(VectorXd, double), VectorXd, double, int);
double fx(VectorXd cj, double x);
double fa(double, double, double(*fx)(VectorXd, double), VectorXd, double, int);
VectorXd FFX(VectorXd);
VectorXd fx(VectorXd);
MatrixXd JacFX(VectorXd);

VectorXd FX(n);
VectorXd aj(n);
VectorXd A0(n);
VectorXd FXold(n);
VectorXd FXnew(n);
VectorXd XX(n);
VectorXd EX(n);
MatrixXd JacobiFX(n,n);
MatrixXd jac(n,n);
VectorXd A1(n);
VectorXd VX(n);

/*
  Expression Of JacobiFX  
MatrixXd
JacFX(VectorXd xx)
{
    int nx;
    int nx1;
    int nx2;
    double tx = 0;
    double ty = 0;

    for (int j = 0; j <= n; j++)
    {
        tx = t1 + j * h2;
        nx = (int)((tx - t1) / h2 + 0.5);
        nx1 = nx - 1;
        nx2 = nx + 1;
        savx = pow(10, tx);
        int z;
        for (int k = 0; k <= n; k++)
        {
            ty = t1 + k * h2;
            savy = pow(10, ty);
            z = (int)((log10(savy) - t1) / h2 + 0.5);
            if (k != j)
            {
                if (k == 0 || k == n)
                {
                    JacobiFX(j, k) = (-0.5) * h2 * trapz(f3, a, b, 100) * (savy - xx(z) * xx(z)) /
                        (savy + xx(z) * xx(z)) / (savy + xx(z) * xx(z));
                }
                else if (k == nx1 || k == nx2)
                {
                    JacobiFX(j, k) = (-1.5) * h2 * trapz(f3, a, b, 100) * (savy - xx(z) * xx(z)) /
                        (savy + xx(z) * xx(z)) / (savy + xx(z) * xx(z));
                }
                else
                {
                    JacobiFX(j, k) = (-h2) * trapz(f3, a, b, 100) * (savy - xx(z) * xx(z)) /
                        (savy + xx(z) * xx(z)) / (savy + xx(z) * xx(z));
                }
            }
            else
            {
                JacobiFX(j, k) = 1.0;
            }
        }
    }
    return JacobiFX;
}
*/


/*  Expression Of JacobiFX  */
MatrixXd
JacFX(VectorXd cj)
{
    aj = cj;
    double s0 = cos(M_PI * 0.5 / n);
    double sn = cos(M_PI * (49 + 0.5) / n);
    double si = 0;
//#pragma omp parallel for num_threads(4)
    for (int i = 0; i < n; i++)
    {
        si = cos(M_PI * (i + 0.5) / n);
        savx = A * aa * pow((A / aa), cos(M_PI * (i + 0.5) / n));
        cout << "jacfx:" << i << endl;
        for (int j = 0; j < n; j++)
        {
            I = j;
            JacobiFX(i,j) = pa(j,si) - ( gausslegend(f5, sn, si, 120) + gausslegend(f5, si, s0, 120) );
        }
    }
 
    return JacobiFX;
}




/*  Expression Of FX  */

VectorXd
FFX(VectorXd cj)
{
    aj = cj;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            XX(i) += cj(j) * cheb(j, cos(M_PI * (i + 0.5) / n));
            
        }
        cout << "XX(i) =   " << XX(i) << endl;
    }
    double si = 0;
    double xi = 0;
    double s0 = cos(M_PI * 0.5 / n);
    double sn = cos(M_PI * (49 + 0.5) / n);

    for (int i = 0; i < n; i++)
    {
        xi = A * aa * pow((A / aa), cos(M_PI * (i + 0.5) / n));
        si = cos(M_PI * (i + 0.5) / n);
        savx = xi;
        EX(i) = gausslegend(f2, sn, si, 120) + gausslegend(f2, si, s0, 120);
        cout << "EX(i) =   " << EX(i) << endl;
    }
    for (int j = 0; j < n; j++)
    {
        FX(j) = XX(j) - EX(j);
        cout << "FX(i) =   "  << FX(j) << endl;
    }
    return FX;
}


/*  Calculate the Partial Differation of cj  */
double
pa(int n, double x)
{
    double d = 0;
    double tj = 0;
    if (n == 0)
    {
        d = 1;
    }
    else
    {
        d = 0;
    }
    tj = cheb(n, x) - 0.5 * d;
    return tj;
}


double
fx(VectorXd cj, double x)
{
    double temp = 0;
    for (int i = 0; i < n; i++)
    {
        temp += cj(i) * cheb(i, x);
    }
    return temp;
}



double
f2(double sy)
{
    Sy = sy;
    return trapz(f3, a, b, 100);
}


double
f3(double t)
{
    return (*intfunc)(savx, t, fx, aj, Sy);
}


double
f(double x, double t, double(*fx)(VectorXd, double), VectorXd cj, double sy)
{
    double y = A * aa * pow((A / aa), sy);
    return (c * y * y * sin(t) * sin(t)) / ((x + y - 2 * sqrt(x * y) * cos(t)) * (1 + c1 * (log(A / (x + y - 2 * sqrt(x * y) * cos(t)))))) * fx(cj,sy) / (y + fx(cj,sy) * fx(cj,sy));
}


double
fa(double x, double t, double(*fx)(VectorXd, double), VectorXd cj, double sy, int i)
{
    double y = A * aa * pow((A / aa), sy);
    return pa(i,sy) * (c * y * y * sin(t) * sin(t)) / ((x + y - 2 * sqrt(x * y) * cos(t)) * (1 + c1 * (log(A / (x + y - 2 * sqrt(x * y) * cos(t)))))) * (y - fx(cj, sy) * fx(cj, sy)) / (y + fx(cj, sy) * fx(cj, sy)) / (y + fx(cj, sy) * fx(cj, sy));
}


double
f4(double t)
{
    return (*intfunc1)(savx, t, fx, aj, Sy, I);
}


double
f5(double sy)
{
    Sy = sy;
    return trapz(f4, a, b, 100);
}




#endif
