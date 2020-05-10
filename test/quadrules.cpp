/*  Realization Methods Used For Solving Numerical Integral  */

#include "quadrules.h"

using namespace std;

/* Chebyshev polynomials:
 * Calculates the <n>-th Chebyshev polynomials at <x>. */
double
cheb(int n, double x)
{
	double a, b, erg;
	if (n == 0)
	{
		erg = 1;
	}
	else if (n == 1)
	{
		erg = x;
	}
	else
	{
		a = 1;
		b = x;
		for (int i = 0; i < n - 1; i++)
		{
			erg = 2.* x * b - a;
			a = b;
			b = erg;
		}
	}
	return erg;
}


/* Legendre polynomial:
 * Calculates the <n>-th Legendre Polynomial at <x>. */
double
evalLeg(int n, double x)
{
    double a, b, erg;
    if (n == 0)
    {
        erg = 1;
    }
    else if (n == 1)
    {
        erg = x;
    }
    else
    {
        a = 1;
        b = x;
        for (int i = 0; i < n - 1; i++)
        {
            double nn = double(i + 2);
            erg = ((2. * nn - 1) * x * b - (nn - 1.) * a) / nn;
            a = b;
            b = erg;
        }
    }
    return erg;
}


/* Gauss-Legendre integration:
 * Gauss-Legendre quadrature of on [-1,1] (http://en.wikipedia.org/wiki/Gaussian_quadrature). The <c_n> weights and nodes are saved in <vnw>.  */
void
GauLeg(int n, double* x)
{
    int nh;
    double z, pnprime, dz;

    nh = (n + 1) / 2;//symmetric nodes

    for (int i = 0; i < nh; i++)//nodes
    {
        z = cos(M_PI * (double(i) + 1. - 0.25) / (double(n) + 0.5)); //initial approximation for n-th node
        dz = 1;
        while (fabs(dz) > 1e-15)//find zero via Newton method
        {
            double pn = evalLeg(n, z);//calculate n-th legendre polynomial
            double pnm1 = evalLeg(n - 1, z);//calculate n-1-th legendre polynomial
            pnprime = double(n) * (z * pn - pnm1) / (z * z - 1.);//derivative of n-th legendre polynomial
            dz = -pn / pnprime;//calculate newton-update
            z += dz;
        }

        x[2 * i] = -z; //store nodes
        x[2 * (n - 1 - i)] = z;
        x[2 * i + 1] = 2.0 / ((1.0 - z * z) * pnprime * pnprime);//calculate and store weights
        x[2 * (n - 1 - i) + 1] = x[2 * i + 1];
    }
}

/* Swap two value <a> and <b>. */
void swap(double& a, double& b) {
    double temp;
    temp = a;
    a = b;
    b = temp;
}

/* Gauss-Jordan algorithm to obtain the solution of <a>.<x>=<b>, where the dimension of <b> is <n>*<m>. <x> replaces <b> and the inverse matrix of <a> replaces <a>. */
void gauss_jordan(double* a, double* b, int n, int m)
{
	int i, icol, irow, j, k, l, ll;
	double big, dum, pivinv;
	int* indxc, * indxr, * ipiv;

	indxc = new int[n];
	indxr = new int[n];
	ipiv = new int[n];

	for (j = 0; j < n; j++) ipiv[j] = 0;

	for (i = 0; i < n; i++) {
		big = 0.0;

		for (j = 0; j < n; j++)
			if (ipiv[j] != 1)
				for (k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if (abs(a[j * n + k]) >= big) {
							big = abs(a[j * n + k]);
							irow = j;
							icol = k;
						}
					}
				}

		++(ipiv[icol]);

		if (irow != icol) {
			for (l = 0; l < n; l++) swap(a[irow * n + l], a[icol * n + l]);
			for (l = 0; l < m; l++) swap(b[irow * m + l], b[icol * m + l]);
		}
		indxr[i] = irow;
		indxc[i] = icol;

		if (a[icol * n + icol] == 0.0) cout << "Matrix inversion: Singular Matrix" << endl;
		pivinv = 1.0 / a[icol * n + icol];
		a[icol * n + icol] = 1.0;
		for (l = 0; l < n; l++) a[icol * n + l] *= pivinv;
		for (l = 0; l < m; l++) b[icol * m + l] *= pivinv;
		for (ll = 0; ll < n; ll++)
			if (ll != icol) {
				dum = a[ll * n + icol];
				a[ll * n + icol] = 0.0;
				for (l = 0; l < n; l++) a[ll * n + l] -= a[icol * n + l] * dum;
				for (l = 0; l < m; l++) b[ll * m + l] -= b[icol * m + l] * dum;
			}
	}

	for (l = n - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l])
			for (k = 0; k < n; k++) {
				swap(a[k * n + indxr[l]], a[k * n + indxc[l]]);
			}
	}

	delete[] indxc;
	delete[] indxr;
	delete[] ipiv;
}


/* Gauss-Jordan inversion of the matrix <a> with size <n>*<x>. Result is stored in <a>. */
void gauss_jordan(double* a, int n) {

	double* b;

	b = new double[n];

	for (int i = 0; i < n; i++) b[i] = 0;

	gauss_jordan(a, b, n, 1);

	delete[] b;

}


/* Print matrix <a> of size <m>*<n> to standard output. */
void printMatrix(double* a, int m, int n) {
	double eps;
	eps = 1e-10;

	for (int i = 0; i < m; i++) {
		cout << "row " << i;
		cout << "( ";
		for (int j = 0; j < n; j++) {
			if (fabs(a[i * m + j]) < 1e-10) a[i * m + j] = 0;
			cout << a[i * m + j] << " ";
		}
		cout << ")" << endl;
	}
	cout << "All entries with an absolute value smaller then " << eps << " have been set to zero in output." << endl;
}


/* Matrix multiplication: a[i,k]*b[k,j]. Size of <a>: <n_row>*<n_sum>; size of <b>: <n_sum>*<n_col>. */
void
multMatrix(double* a, double* b, double* erg, int n_row, int n_sum, int n_col) {

	for (int i = 0; i < n_row; i++) {//rows
		for (int j = 0; j < n_col; j++) {//columns
			erg[i * n_col + j] = 0;
			for (int k = 0; k < n_sum; k++) erg[i * n_col + j] += a[i * n_sum + k] * b[k * n_col + j];
		}
	}
}
