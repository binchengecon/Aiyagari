#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace ::std;
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
// rowa and cola are no of rows and columns
// of matrix A
// rowb and colb are no of rows and columns
// of matrix B

const int cola = 2, rowa = 2, colb = 2, rowb = 2;

// Function to computes the Kronecker Product
// of two matrices
void Kroneckerproduct(int A[rowa][cola], int B[rowb][colb], int C[rowa * rowb][cola * colb])
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
                    std::cout << C[k + i * rowb][j * colb + l] << " ";
                }
            }
            std::cout << "\n";
        }
    }
}

int main()
{

    // Driver Code
    int A[2][2] = {{1, 2}, {3, 4}},
        B[2][2] = {{0, 5}, {6, 7}};
    int C[4][4];
    int D[4][4];

    Kroneckerproduct(A, B, C);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << C[i][j] << " ";
        }
        std::cout << "\n";
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << C[i][j] << " ";
        }
        std::cout << "\n";
    }

        return 0;
}
