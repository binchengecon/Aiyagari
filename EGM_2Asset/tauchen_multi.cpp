// Contributor: Bin (2022)

// Purpose:
// Approximate a multivariate AR(1) process with a N^dim-states (node) first order process

// Detail:
// Suppose we have a multivariate AR(1) as Z_t = A1 +A2 Z_{t-1} + eps_t, eps_t~N(0,Sigma)
// Another form would be the following set of equations:
// r_t = (1-rho1)*mu_r + rho1* r_{t-1}+eps1_t
// y_t = (1-rho2)*mu_y + rho2* y_{t-1}+eps2_t
// Sigma = [ sigma1^2   corr*sigma1*sigma2
//  		corr*sigma1*sigma2		sigma2^2]
// Then we utilize Vec(Var) = (I_{n^2}-A otimes A)^{-1}Vec(Sigma)
// To approximate Z_t, one sets up a finite grid S consisting of Q=n1*n2 possible states.
//  R={r_0,r_1,..,r_{n1-1}}
//  Y={y_0,y_1,..,y_{n2-1}}
// 1. sigma_r = sqrt(Vec(Var)[0])
// 2. sigma_y = sqrt(Vec(Var)[3])
// 3. r_0 = mu_r - m*sigma_r, r_{n1-1}=mu_r+m*sigma_r, delta_r = (r_{n1-1}-r_0)/(n1-1)
// 4. y_0 = mu_y - m*sigma_y, y_{n2-1}=mu_y+m*sigma_y, delta_y = (y_{n2-1}-y_0)/(n2-1)

// S = {s1,s2,...,sQ} where s_k = (r_{k/n1},y_{k mod n1})
// Transition matrix P would be
// P_{i,j} = P(M_{t+1}=sj|M_t=si)
// And approximation is that
// P_{i,j} = P(Z_{t+1} \in Vj|Z_t=si) s.t. si \in V_i and U_{i=1}^{Q}Vi = R^m
// 							jr=0,jy=0				(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&					(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
// 							jr=0,0<jy<n2-1			(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
// 							jr=0,jy=n2-1			(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon
// P_{(ir,iy),(jr,jy)}={    0<jr<n1-1,jy=0          r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2						 }
//							0<jr<n1-1,0<jy<n2-1		r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
//							0<jr<n1-1,jy=n2-1		r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon
//							jr=n1-1,jy=0			r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
//							jr=n1-1,0<jy<n2-1		r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
// 							jr=n1-1,jy=n2-1			r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon

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

double CDFSTDNormal(double x)
{
	// Function CDFSTDNormal: computes the standard normal CDF using Abramowiz and Stegun (1964) approximation

	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x) / sqrt(2.0);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

	return 0.5 * (1.0 + sign * y);
}

double CDFSTDNormal_2D_r1y1(double ah, double ak, double r)

//****************************************************************************80
//
//  Purpose:
//
//    BIVNOR computes the bivariate normal CDF.
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

double CDFSTDNormal_2D_r1y0(double ah, double ak, double r)
{
	return (1 - CDFSTDNormal(ah)) - CDFSTDNormal_2D_r1y1(ah, ak, r);
}

double CDFSTDNormal_2D_r0y1(double ah, double ak, double r)
{
	return (1 - CDFSTDNormal(ak)) - CDFSTDNormal_2D_r1y1(ah, ak, r);
}

double CDFSTDNormal_2D_r0y0(double ah, double ak, double r)
{
	return 1 - CDFSTDNormal_2D_r1y1(ah, ak, r) - CDFSTDNormal_2D_r0y1(ah, ak, r) - CDFSTDNormal_2D_r1y1(ah, ak, r);
}

const int cola = 2, rowa = 2, colb = 2, rowb = 2;

// Function to computes the Kronecker Product
// of two matrices
void Kroneckerproduct(double A[rowa][cola], double B[rowb][colb], double C[rowa * rowb][cola * colb])
{

	// i loops till rowa
	for (int i = 0; i < rowa; i++)
	{

		// k loops till rowb
		for (int j = 0; j < cola; j++)
		{

			// j loops till cola
			for (int k = 0; k < rowb; k++)
			{

				// l loops till colb
				for (int l = 0; l < colb; l++)
				{

					// Each element of matrix A is
					// multiplied by whole Matrix B
					// resp and stored as Matrix C
					C[k + i * rowb][j * colb + l] = A[i][j] * B[k][l];
					// std::cout << C[k + i * rowb][j * colb + l] << " ";
				}
			}
			// std::cout << "\n";
		}
	}
}

// void MatProduct(double M[4][4], double MI[4][4],double I[4][4])
void MatProductMat(double M[4][4], double MI[4][4])
{

	int i, j, k;
	double temp;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			temp = 0;
			for (k = 0; k < 4; k++)
			{
				temp += M[i][k] * MI[k][j];
			}
			std::cout << temp << " ";
			// I[i][j]=temp;
			temp = 0;
		}
		std::cout << "\n";
	}
}

void eye(double I[4][4])
{
	int i, j;
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			if (i == j)
			{
				I[i][j] = 1;
			}
			else
			{
				I[i][j] = 0;
			}
		}
	}
}
// void MatProduct(double M[4][4], double MI[4][4],double I[4][4])
void MatDeductMat(double M[4][4], double MI[4][4], double MII[4][4])
{

	int i, j, k;
	double temp;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			MII[i][j] = M[i][j] - MI[i][j];
		}
	}
}

// void MatProduct(double M[4][4], double MI[4][4],double I[4][4])
void MatProductVec(double M[4][4], double VR[4], double VL[4])
{

	int i, j, k;
	double temp;

	for (i = 0; i < 4; i++)
	{

		temp = 0;
		for (k = 0; k < 4; k++)
		{
			temp += M[i][k] * VR[k];
		}
		// std::cout << temp << " ";
		VL[i] = temp;
		temp = 0;

		// std::cout << "\n";
	}
}

bool InverseMat44(double m[16], double invOut[16])
{
	double inv[16], det;
	int i;

	inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
	inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
	inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
	inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
	inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
	inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
	inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
	inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
	inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
	inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
	inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
	inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
	inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
	inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
	inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
	inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
	if (det == 0)
		return false;
	det = 1.0 / det;

	for (i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	return true;
}

void Vec(double M[4][4], double M_vec[16])
{
	int i, j;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			M_vec[i * 4 + j] = M[i][j];
		}
	}
}

void Vec2(double M[2][2], double M_vec[4])
{
	int i, j;

	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
		{
			M_vec[i * 2 + j] = M[i][j];
		}
	}
}

void deVec(double M_vec[16], double M[4][4])
{
	int i, j;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			M[i][j] = M_vec[i * 4 + j];
		}
	}
}

template <size_t Nodes>
void tauchenfun2D(double rho1, double rho2, double m, double amu1, double amu2, double sigma1, double sigma2, double corr, double grid1[Nodes * Nodes], double gird2[Nodes * Nodes], double Ptransi[Nodes * Nodes][Nodes * Nodes])
{

	int i, j, k;

	double AR[2][2] = {{rho1, 0}, {0, rho2}},
		   Sigma[2][2] = {{sigma1 * sigma1, corr * sigma1 * sigma2}, {corr * sigma1 * sigma2, sigma2 * sigma2}}, I[4][4];
	double ARAR[4][4], IDARAR[4][4], IDARAR_inv[4][4], Sigma_vec[4], Var_Vec[4];
	double IDARAR_Vec[16], IDARAR_Vec_inv[16];

	Kroneckerproduct(AR, AR, ARAR);
	eye(I);
	MatDeductMat(I, ARAR, IDARAR);
	Vec(IDARAR, IDARAR_Vec);
	InverseMat44(IDARAR_Vec, IDARAR_Vec_inv);
	deVec(IDARAR_Vec_inv, IDARAR_inv);
	Vec2(Sigma, Sigma_vec);
	MatProductVec(IDARAR_inv, Sigma_vec, Var_Vec);

	// First lets compute unconditiontal variance of yt

	double var1t = (Var_Vec[0] * Var_Vec[0]) / (1 - rho1 * rho1);
	double var2t = (Var_Vec[3] * Var_Vec[3]) / (1 - rho2 * rho2);

	// Compute stddev of yt

	double std1t = sqrt(var1t);
	double std2t = sqrt(var2t);

	// std::cout << "value of std_e="<<stdyt <<"\n";
	// Define maximum and minimum grid point

	double y1nodes[Nodes], y2nodes[Nodes];

	y1nodes[Nodes - 1] = m * std1t;
	y1nodes[0] = -y1nodes[Nodes - 1];

	// Define interior nodes

	double y1nodesinterval = (y1nodes[Nodes - 1] - y1nodes[0]) / ((Nodes - 1) * 1.0);

	for (i = 1; i < (Nodes - 1); i++)
	{
		y1nodes[i] = y1nodes[i - 1] + y1nodesinterval;
	}

	for (i = 0; i < Nodes; i++)
	{
		y1nodes[i] = y1nodes[i] + amu1;
	}

	y2nodes[Nodes - 1] = m * std1t;
	y2nodes[0] = -y2nodes[Nodes - 1];

	// Define interior nodes

	double y2nodesinterval = (y2nodes[Nodes - 1] - y2nodes[0]) / ((Nodes - 1) * 1.0);

	for (i = 1; i < (Nodes - 1); i++)
	{
		y2nodes[i] = y2nodes[i - 1] + y2nodesinterval;
	}

	for (i = 0; i < Nodes; i++)
	{
		y2nodes[i] = y2nodes[i] + amu2;
	}

	for (i = 0; i < Nodes; i++)
	{
		for (j = 0; j < Nodes; j++)
		{
			grid1[i * Nodes + j] = y1nodes[i];
			grid2[i * Nodes + j] = y2nodes[j];
		}
	}

	// Computing transition probability matrix

	double transitionMat[Nodes * Nodes][Nodes * Nodes], temp, tempr, tempy;

	int ir, iy, jr, jy, index_ry, index_ry_next;
	for (ir = 0; ir < Nodes; ir++)
	{
		for (iy = 0; iy < Nodes; iy++)
		{

			// starting at state = ir, iy
			tempr = (1 - rho1) * amu1 + rho1 * y1nodes[ir];
			tempy = (1 - rho2) * amu2 + rho2 * y2nodes[iy];
			index_ry = ir * Nodes + iy;

			// 							jr=0,jy=0				(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&					(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
			// 							jr=0,0<jy<n2-1			(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
			// 							jr=0,jy=n2-1			(1-rho1)*mu_r + rho1*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon

			jr = 0;

			for (jy = 1; jy < (Nodes - 1); jy++)
			{
				temp = 0;
				index_ry_next = jr * Nodes + jy;
				temp = CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t, corr);
				temp -= CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t, corr);
				transitionMat[index_ry][index_ry_next] = temp;
				temp = 0;
				//
			}

			jy = 0;

			temp = 0;
			index_ry_next = jr * Nodes + jy;
			temp = CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t, corr);
			transitionMat[index_ry][index_ry_next] = temp;
			temp = 0;

			jy = (Nodes - 1);

			temp = 0;
			index_ry_next = jr * Nodes + jy;
			temp = CDFSTDNormal_2D_r0y1((y1nodes[jr] + y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t, corr);
			transitionMat[index_ry][index_ry_next] = temp;
			temp = 0;

			// P_{(ir,iy),(jr,jy)}={    0<jr<n1-1,jy=0          r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2						 }
			//							0<jr<n1-1,0<jy<n2-1		r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
			//							0<jr<n1-1,jy=n2-1		r_j-delta_r/2 < (1-rho2)*mu_r + rho2*r_i+epsilon <  r_j+delta_r/2 &&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon

			for (jr = 1; jr < (Nodes - 1); jr++)
			{

				for (jy = 1; jy < (Nodes - 1); jy++)
				{
					temp = 0;
					index_ry_next = jr * Nodes + jy;
					temp = CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t, corr);
					temp -= CDFSTDNormal_2D_r0y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t, corr);
					transitionMat[index_ry][index_ry_next] = temp;
					temp = 0;
					//
				}

				jy = 0;

				temp = 0;
				index_ry_next = jr * Nodes + jy;
				temp = CDFSTDNormal_2D_r0y0((y1nodes[jr] + y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t, corr);
				temp -= CDFSTDNormal_2D_r0y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t, corr);
				transitionMat[index_ry][index_ry_next] = temp;
				temp = 0;

				jy = (Nodes - 1);

				temp = 0;
				index_ry_next = jr * Nodes + jy;
				temp = CDFSTDNormal_2D_r0y1((y1nodes[jr] + y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t, corr);
				temp -= CDFSTDNormal_2D_r0y1((y1nodes[jr] - y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t, corr);
				transitionMat[index_ry][index_ry_next] = temp;
				temp = 0;
			}

			//							jr=n1-1,jy=0			r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	(1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
			//							jr=n1-1,0<jy<n2-1		r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon <  y_j+delta_y/2
			// 							jr=n1-1,jy=n2-1			r_j-delta_r/2 <	(1-rho1)*mu_r + rho1*r_i+epsilon	&&	y_j-delta_y/2 < (1-rho2)*mu_y + rho2*y_i+epsilon

			jr = Node - 1;

			for (jy = 1; jy < (Nodes - 1); jy++)
			{
			temp = 0;
			index_ry_next = jr * Nodes + jy;
			temp = CDFSTDNormal_2D_r1y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] + y2nodesinterval / 2 - tempy) / std2t, corr);
			temp -= CDFSTDNormal_2D_r1y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t, corr);
			transitionMat[index_ry][index_ry_next] = temp;
			temp = 0;				//
			}

			jy = 0;

			temp = 0;
			index_ry_next = jr * Nodes + jy;
			temp = CDFSTDNormal_2D_r1y0((y1nodes[jr] - y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t, corr);
			transitionMat[index_ry][index_ry_next] = temp;
			temp = 0;

			jy = (Nodes - 1);

			temp = 0;
			index_ry_next = jr * Nodes + jy;
			temp -= CDFSTDNormal_2D_r1y1((y1nodes[jr] - y1nodesinterval / 2 - tempr) / std1t, (y2nodes[jy] - y2nodesinterval / 2 - tempy) / std2t, corr);
			transitionMat[index_ry][index_ry_next] = temp;
			temp = 0;
		}

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