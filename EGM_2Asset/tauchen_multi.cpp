// Contributor: Bin (2022)

// Purpose:
// Approximate a multivariate AR(1) process with a N^dim-states (node) first order process

// Detail:
// Suppose we have amultivariate AR(1) as Z_t = A1 +A2 Z_{t-1} + eps_t, eps_t~N(0,Sigma)
// To approximate Z_t, one sets up a finite grid S consisting of Q possible states.
//  S = {s1,s2,...,sQ} \in R^m
// Transition matrix P would be
// P_{i,j} = P(M_{t+1}=sj|M_t=si)
// And approximation is that
// P_{i,j} = P(Z_{t+1} \in Vj|Z_t=si) s.t. si \in V_i and U_{i=1}^{Q}Vi = R^m

// AR(1) as y_t = (1-A[2])*mu_y + A[2]*y_{t-1}+ epsilon_t, where epsilon_t is N(0,sigma^2)
// We use n nodes to approxiamate this continuously valued markov chain
// The space Y={y_0,y_2,..,y_{n-1}}
// The transition matrix is P={p_{i,j}}_{0<= i,j <n}
// 1. sigma_y = sigma/sqrt(1-A[2]^2)
// 2. y_0 = mu_y - m*sigma_y, y_{n-1}=mu_y+m*sigma_y, delta_y = (y_{n-1}-y_0)/(n-1)
// 3.                               0=j,                    (1-A[2])*mu_y + A[2]*y_i+epsilon <  y_j+delta_y/2
//              p_{i,j} = {         0<j<n-1, y_j-delta_y/2 < (1-A[2])*mu_y + A[2]*y_i+epsilon <  y_j+delta_y/2 }
//                                  n-1=j,                  y_j-delta_y/2 < (1-A[2])*mu_y + A[2]*y_i+epsilon

//--------Functions definition:
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace std;

double gauss(double t)

//****************************************************************************80
//
//  Purpose:
//
//    GAUSS returns the area of the lower tail of the normal curve.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the evaluation point.
//
//    Output, double GAUSS, the lower normal tail area.
//
{
	double value;

	value = (1.0 + erf(t / sqrt(2.0))) / 2.0;

	return value;
}

double CDFSTDNormal_2D(double ah, double ak, double r)

//****************************************************************************80
//
//  Purpose:
//
//    BIVNOR computes the bivariate normal CDF.
//
//  Discussion:
//
//    BIVNOR computes the probability for two normal variates X and Y
//    whose correlation is R, that AH <= X and AK <= Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    Original FORTRAN77 version by Thomas Donnelly.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Thomas Donnelly,
//    Algorithm 462: Bivariate Normal Distribution,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 10, page 638.
//
//  Parameters:
//
//    Input, double AH, AK, the lower limits of integration.
//
//    Input, double R, the correlation between X and Y.
//
//    Output, double BIVNOR, the bivariate normal CDF.
//
//  Local Parameters:
//
//    Local, int IDIG, the number of significant digits
//    to the right of the decimal point desired in the answer.
//
{
	double a2;
	double ap;
	double b;
	double cn;
	double con;
	double conex;
	double ex;
	double g2;
	double gh;
	double gk;
	double gw;
	double h2;
	double h4;
	int i;
	static int idig = 15;
	int is;
	double rr;
	double s1;
	double s2;
	double sgn;
	double sn;
	double sp;
	double sqr;
	double t;
	static double twopi = 6.283185307179587;
	double w2;
	double wh;
	double wk;

	b = 0.0;

	gh = gauss(-ah) / 2.0;
	gk = gauss(-ak) / 2.0;

	if (r == 0.0)
	{
		b = 4.00 * gh * gk;
		b = fmax(b, 0.0);
		b = fmin(b, 1.0);
		return b;
	}

	rr = (1.0 + r) * (1.0 - r);

	if (rr < 0.0)
	{
		cerr << "\n";
		cerr << "BIVNOR - Fatal error!\n";
		cerr << "  1 < |R|.\n";
		exit(0);
	}

	if (rr == 0.0)
	{
		if (r < 0.0)
		{
			if (ah + ak < 0.0)
			{
				b = 2.0 * (gh + gk) - 1.0;
			}
		}
		else
		{
			if (ah - ak < 0.0)
			{
				b = 2.0 * gk;
			}
			else
			{
				b = 2.0 * gh;
			}
		}
		b = fmax(b, 0.0);
		b = fmin(b, 1.0);
		return b;
	}

	sqr = sqrt(rr);

	if (idig == 15)
	{
		con = twopi * 1.0E-15 / 2.0;
	}
	else
	{
		con = twopi / 2.0;
		for (i = 1; i <= idig; i++)
		{
			con = con / 10.0;
		}
	}
	//
	//  (0,0)
	//
	if (ah == 0.0 && ak == 0.0)
	{
		b = 0.25 + asin(r) / twopi;
		b = fmax(b, 0.0);
		b = fmin(b, 1.0);
		return b;
	}
	//
	//  (0,nonzero)
	//
	if (ah == 0.0 && ak != 0.0)
	{
		b = gk;
		wh = -ak;
		wk = (ah / ak - r) / sqr;
		gw = 2.0 * gk;
		is = 1;
	}
	//
	//  (nonzero,0)
	//
	else if (ah != 0.0 && ak == 0.0)
	{
		b = gh;
		wh = -ah;
		wk = (ak / ah - r) / sqr;
		gw = 2.0 * gh;
		is = -1;
	}
	//
	//  (nonzero,nonzero)
	//
	else if (ah != 0.0 && ak != 0.0)
	{
		b = gh + gk;
		if (ah * ak < 0.0)
		{
			b = b - 0.5;
		}
		wh = -ah;
		wk = (ak / ah - r) / sqr;
		gw = 2.0 * gh;
		is = -1;
	}

	for (;;)
	{
		sgn = -1.0;
		t = 0.0;

		if (wk != 0.0)
		{
			if (fabs(wk) == 1.0)
			{
				t = wk * gw * (1.0 - gw) / 2.0;
				b = b + sgn * t;
			}
			else
			{
				if (1.0 < fabs(wk))
				{
					sgn = -sgn;
					wh = wh * wk;
					g2 = gauss(wh);
					wk = 1.0 / wk;

					if (wk < 0.0)
					{
						b = b + 0.5;
					}
					b = b - (gw + g2) / 2.0 + gw * g2;
				}
				h2 = wh * wh;
				a2 = wk * wk;
				h4 = h2 / 2.0;
				ex = exp(-h4);
				w2 = h4 * ex;
				ap = 1.0;
				s2 = ap - ex;
				sp = ap;
				s1 = 0.0;
				sn = s1;
				conex = fabs(con / wk);

				for (;;)
				{
					cn = ap * s2 / (sn + sp);
					s1 = s1 + cn;

					if (fabs(cn) <= conex)
					{
						break;
					}
					sn = sp;
					sp = sp + 1.0;
					s2 = s2 - w2;
					w2 = w2 * h4 / sp;
					ap = -ap * a2;
				}
				t = (atan(wk) - wk * s1) / twopi;
				b = b + sgn * t;
			}
		}
		if (0 <= is)
		{
			break;
		}
		if (ak == 0.0)
		{
			break;
		}
		wh = -ak;
		wk = (ah / ak - r) / sqr;
		gw = 2.0 * gk;
		is = 1;
	}

	b = fmax(b, 0.0);
	b = fmin(b, 1.0);

	return b;
}

const int cola = 2, rowa = 2, colb = 2, rowb = 2;

// Function to computes the Kronecker Product
// of two matrices
void Kroneckerproduct(int A[rowa][cola], int B[rowb][colb], int C[rowa * rowb][cola * colb])
{

	// i loops till rowa
	for (int i = 0; i < rowa; i++)
	{

		// k loops till rowb
		for (int k = 0; k < rowb; k++)
		{

			// j loops till cola
			for (int j = 0; j < cola; j++)
			{

				// l loops till colb
				for (int l = 0; l < colb; l++)
				{

					// Each element of matrix A is
					// multiplied by whole Matrix B
					// resp and stored as Matrix C
					C[i + l + 1][j + k + 1] = A[i][j] * B[k][l];
					cout << C[i + l + 1][j + k + 1] << " ";
				}
			}
			std::cout << endl;
		}
	}
}

template <size_t Nodes>
void tauchenfun2D(double A[2], double m, double amu[2], double sigma[2][2], double grid[Nodes], double Ptransi[Nodes][Nodes])

{

	int i, j, k;

	// First lets compute unconditiontal variance of yt

	double varyt = (sigma * sigma) / (1 - A[2] * A[2]);

	// Compute stddev of yt

	double stdyt = sqrt(varyt);
	// std::cout << "value of std_e="<<stdyt <<"\n";
	// Define maximum and minimum grid point

	double ynodes[Nodes];
	ynodes[Nodes - 1] = m * stdyt;
	ynodes[0] = -ynodes[Nodes - 1];

	// Define interior nodes

	double ynodesinterval = (ynodes[Nodes - 1] - ynodes[0]) / ((Nodes - 1) * 1.0);

	for (i = 1; i < (Nodes - 1); i++)
	{
		ynodes[i] = ynodes[i - 1] + ynodesinterval;
	}

	for (i = 0; i < Nodes; i++)
	{
		ynodes[i] = ynodes[i] + amu;
	}

	for (i = 0; i < Nodes; i++)
	{
		grid[i] = ynodes[i];
	}

	// Computing transition probability matrix

	double transitionMat[Nodes][Nodes];

	for (j = 0; j < Nodes; j++)
	{
		for (k = 1; k < (Nodes - 1); k++)
		{

			transitionMat[j][k] = CDFSTDNormal((ynodes[k] - (1 - A[2]) * amu - A[2] * ynodes[j] + ynodesinterval / 2.0) / sigma) - CDFSTDNormal((ynodes[k] - (1 - A[2]) * amu - A[2] * ynodes[j] - ynodesinterval / 2.0) / sigma);
		}

		transitionMat[j][0] = CDFSTDNormal((ynodes[0] - (1 - A[2]) * amu - A[2] * ynodes[j] + ynodesinterval / 2.0) / sigma);
		transitionMat[j][Nodes - 1] = 1.0 - CDFSTDNormal((ynodes[Nodes - 1] - (1 - A[2]) * amu - A[2] * ynodes[j] - ynodesinterval / 2.0) / sigma);
	}

	for (j = 0; j < Nodes; j++)
	{
		for (k = 0; k < (Nodes); k++)
		{
			Ptransi[j][k] = transitionMat[j][k];
		}
	}
}

// Gaussian Quadratures.
typedef double DP;
void gauher(double x[], double w[], int ngrid)
{
	const DP EPS = 1.0e-14, PIM4 = 0.7511255444649425;
	const int MAXIT = 10;
	int i, its, j, m;
	DP p1 = 0, p2 = 0, p3 = 0, pp = 0, z = 0, z1 = 0;

	int n = ngrid;
	m = (n + 1) / 2;
	for (i = 0; i < m; i++)
	{
		if (i == 0)
		{
			z = sqrt(DP(2 * n + 1)) - 1.85575 * pow(DP(2 * n + 1), -0.16667);
		}
		else if (i == 1)
		{
			z -= 1.14 * pow(DP(n), 0.426) / z;
		}
		else if (i == 2)
		{
			z = 1.86 * z - 0.86 * x[0];
		}
		else if (i == 3)
		{
			z = 1.91 * z - 0.91 * x[1];
		}
		else
		{
			z = 2.0 * z - x[i - 2];
		}
		for (its = 0; its < MAXIT; its++)
		{
			p1 = PIM4;
			p2 = 0.0;
			for (j = 0; j < n; j++)
			{
				p3 = p2;
				p2 = p1;
				p1 = z * sqrt(2.0 / (j + 1)) * p2 - sqrt(DP(j) / (j + 1)) * p3;
			}
			pp = sqrt(DP(2 * n)) * p2;
			z1 = z;
			z = z1 - p1 / pp;
			if (fabs(z - z1) <= EPS)
				break;
		}
		if (its >= MAXIT)
			printf("too many iterations in gauher");
		x[i] = z;
		x[n - 1 - i] = -z;
		w[i] = 2.0 / (pp * pp);
		w[n - 1 - i] = w[i];
	}
}