#define _USE_MATH_DEFINES 
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <ctime>
#include <vector>

using namespace std;
using namespace Eigen;

/*  Initialize The Parameters  */
const int n = 100;
const double a = 0;
const double b = M_PI;
const double aa = 0.1;
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
double f(double , double , double );
double f2(double );
double f3(double );
double savx;
double savy;
double (*intfunc)(double, double, double);
void FFX(VectorXd );
void JacFX(VectorXd );
VectorXd FX(n + 1);
VectorXd XX(n + 1);
MatrixXd JacobiFX(n + 1, n + 1);
VectorXd X1(n + 1);
VectorXd VX(n + 1);

/*  Quadrules Of Functions */
template <class T>
double trapz(T& func, const double a, const double b)
{
    double h = (b - a) / n;
    double s = 0;
    double sum = 0;
    for (int i = 1; i < n; i++)  s += func(a + i * h);
    sum += h * (func(a) + 2 * s + func(b)) / 2;
    return sum;
}

template <class T>
double trapz(T& func, const double xi, const double a, const double b)
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



int main()
{
    //clock_t start = clock();
    intfunc = f;
   
    for (int i = 0; i <= n; i++)
    {
        XX(i) = 110.0;
    }
   
    
    
    for (int j = 0; j <= 1000; j++)
    {
        FFX(XX);
        JacFX(XX);
        VX = JacobiFX.inverse()  * FX;
        if (VX.norm() < 1e-8)
            break;
        X1 = XX - VX;
        XX = X1;
        cout << VX.norm() << endl;
    }
    

    ofstream out("X1.dat");
    for (int i = 0; i <= n; i++)
    {
        out << X1(i) << "\t" << endl;   
    }
    out.close();

    // clock_t finish = clock();
    // double TotalTime = (double)(finish - start) / CLOCKS_PER_SEC;
    // cout << "\n程序的运行时间为" << TotalTime << "秒！" << endl;

    return 0;
}


/*  Expression Of FX  */
void
FFX(VectorXd XX)
{
    double tx = 0;
    for (int j = 0; j <= n; j++)
    {
        tx = t1 + j * h2;
        savx = pow(10, tx);
        FX(j) = XX(j) - trapz(f2, tx, t1, t2);
    }
    
}


/*  Expression Of JacobiFX  */
void
JacFX(VectorXd XX)
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
                    JacobiFX(j, k) = (-0.5) * h2 * trapz(f3, a, b) * (savy - XX(z) * XX(z)) /
                        (savy + XX(z) * XX(z)) / (savy + XX(z) * XX(z));
                }
                else if (k == nx1 || k == nx2)
                {
                    JacobiFX(j, k) = (-1.5) * h2 * trapz(f3, a, b) * (savy - XX(z) * XX(z)) /
                        (savy + XX(z) * XX(z)) / (savy + XX(z) * XX(z));
                }
                else
                {
                    JacobiFX(j, k) = (-h2) * trapz(f3, a, b) * (savy - XX(z) * XX(z)) /
                        (savy + XX(z) * XX(z)) / (savy + XX(z) * XX(z));
                }
            }
            else
            {
                JacobiFX(j, k) = 1.0;
            }
        }
    }
    
}



double 
f2(double y)
{
    savy = pow(10, y);
    return trapz(f3, a, b);
}


double 
f3(double t)
{
    return (*intfunc)(savx, savy, t);
}


double 
f(double x, double y, double t)
{
    int z = (int)((log10(y) - t1) / h2 + 0.5);
    return (c * y * y * sin(t) * sin(t)) / ((x + y - 2 * sqrt(x * y) * cos(t)) * (1 + c1 * (log(A / (x + y - 2 * sqrt(x * y) * cos(t)))))) * XX(z) / (y + XX(z) * XX(z));
}




/*
double
funcT(double t)
{
    ft = (c * y * y * sin(t) * sin(t)) / ((x + y - 2 * sqrt(x * y) * cos(t)) * (1 + c1 * (log(A / (x + y - 2 * sqrt(x * y) * cos(t))))) );
    return ft;
}


double
funcY(double y)
{
    int z = (int)( (y - aa) / h1 + 0.5);
    fy = trapz(funcT, a, b, n) * XX(z) / (y + XX(z) * XX(z));
    return  fy;
}
*/



