// Update on the whole algorithm: next period and current period return rate timing is different from what I thought previously.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <sys/time.h>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <direct.h>
#include <iomanip>
#include <cmath>
#include <limits>
#include <windows.h>
#include <string>
#include <sstream>

const int size_k = 500; // number of grid points
const int size_y = 7;   // number of productivity classes
const int size_j = 200; // 1000, 4000 (stuck at 2.2e-004) not working

#define ifulldim (size_k * size_y * size_y * (size_j + 1))
#define iydim (size_k * size_y * size_y)

#define inx(igridindex, jclassindex, jjclassindex, kclassindex) (((kclassindex) * (size_k * size_y * size_y)) + ((jjclassindex) * (size_k * size_y)) + ((jclassindex) * (size_k)) + (igridindex))
#define inx2(igridindex, jclassindex, jjclassindex) (((jjclassindex) * (size_k * size_y)) + ((jclassindex) * (size_k)) + (igridindex))

const double kmin = 0.0;
const double kmax = 500.0;

const double betapar = 0.9;
const double alphapar = 0.36;
const double deltapar = 0.08;
const double rhopar = 3.0;
const double labor = 1.0219882;

const double epsV = 1.0e-8;
const double epsdist = 1.0e-7;
const double epsK = 1.0e-6;
const double relaxsK = 0.005;
const double relaxVF = 0.995;

// grid constants
const double scale1 = 1.6;
const double grmin = (kmin / scale1) - 1.0;
const double exponen = log((kmax / scale1) - grmin) / (size_k - 1);

// const double risk_states[7] = {-0.030619, -0.020412, -0.010206, 0.000000, 0.010206, 0.020412, 0.030619};
// const double risk_trans[7][7] = {
//     {0.026240, 0.152924, 0.361483, 0.328567, 0.114742, 0.015266, 0.000778},
//     {0.016044, 0.114742, 0.328567, 0.361483, 0.152924, 0.024700, 0.001539},
//     {0.009452, 0.082835, 0.287445, 0.382789, 0.196114, 0.038437, 0.002929},
//     {0.005362, 0.057531, 0.242024, 0.390166, 0.242024, 0.057531, 0.005362},
//     {0.002929, 0.038437, 0.196114, 0.382789, 0.287445, 0.082835, 0.009452},
//     {0.001539, 0.024700, 0.152924, 0.361483, 0.328567, 0.114742, 0.016044},
//     {0.000778, 0.015266, 0.114742, 0.328567, 0.361483, 0.152924, 0.026240}};

// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.006262, 0.060934, 0.242398, 0.382924, 0.241063, 0.060262, 0.006157},
//     {0.006245, 0.060822, 0.242175, 0.382924, 0.241285, 0.060374, 0.006175},
//     {0.006227, 0.060710, 0.241953, 0.382925, 0.241508, 0.060486, 0.006192},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006192, 0.060486, 0.241508, 0.382925, 0.241953, 0.060710, 0.006227},
//     {0.006175, 0.060374, 0.241285, 0.382924, 0.242175, 0.060822, 0.006245},
//     {0.006157, 0.060262, 0.241063, 0.382924, 0.242398, 0.060934, 0.006262}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.006262, 0.060934, 0.242398, 0.382924, 0.241063, 0.060262, 0.006157},
//     {0.006245, 0.060822, 0.242175, 0.382924, 0.241285, 0.060374, 0.006175},
//     {0.006227, 0.060710, 0.241953, 0.382925, 0.241508, 0.060486, 0.006192},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006192, 0.060486, 0.241508, 0.382925, 0.241953, 0.060710, 0.006227},
//     {0.006175, 0.060374, 0.241285, 0.382924, 0.242175, 0.060822, 0.006245},
//     {0.006157, 0.060262, 0.241063, 0.382924, 0.242398, 0.060934, 0.006262}};

// const double p_e = 0.000001;
// const double std_e = 0.01;
// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241730, 0.382925, 0.241730, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241730, 0.382925, 0.241730, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210}};

// const double p_e = 0.0000001;
// const double std_e = 0.01;

// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.0062096854, 0.0605975845, 0.2417303731, 0.3829249254, 0.2417302395, 0.0605975173, 0.0062096749},
//     {0.0062096836, 0.0605975733, 0.2417303508, 0.3829249254, 0.2417302618, 0.0605975285, 0.0062096766},
//     {0.0062096819, 0.0605975621, 0.2417303286, 0.3829249254, 0.2417302840, 0.0605975397, 0.0062096784},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096784, 0.0605975397, 0.2417302840, 0.3829249254, 0.2417303286, 0.0605975621, 0.0062096819},
//     {0.0062096766, 0.0605975285, 0.2417302618, 0.3829249254, 0.2417303508, 0.0605975733, 0.0062096836},
//     {0.0062096749, 0.0605975173, 0.2417302395, 0.3829249254, 0.2417303731, 0.0605975845, 0.0062096854}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.0062096854, 0.0605975845, 0.2417303731, 0.3829249254, 0.2417302395, 0.0605975173, 0.0062096749},
//     {0.0062096836, 0.0605975733, 0.2417303508, 0.3829249254, 0.2417302618, 0.0605975285, 0.0062096766},
//     {0.0062096819, 0.0605975621, 0.2417303286, 0.3829249254, 0.2417302840, 0.0605975397, 0.0062096784},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096784, 0.0605975397, 0.2417302840, 0.3829249254, 0.2417303286, 0.0605975621, 0.0062096819},
//     {0.0062096766, 0.0605975285, 0.2417302618, 0.3829249254, 0.2417303508, 0.0605975733, 0.0062096836},
//     {0.0062096749, 0.0605975173, 0.2417302395, 0.3829249254, 0.2417303731, 0.0605975845, 0.0062096854}};

// const double p_e = 0.00000001;
// const double std_e = 0.01;

// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.0062096806, 0.0605975542, 0.2417303130, 0.3829249254, 0.2417302996, 0.0605975475, 0.0062096796},
//     {0.0062096805, 0.0605975531, 0.2417303108, 0.3829249254, 0.2417303018, 0.0605975486, 0.0062096798},
//     {0.0062096803, 0.0605975520, 0.2417303085, 0.3829249254, 0.2417303041, 0.0605975497, 0.0062096799},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096799, 0.0605975497, 0.2417303041, 0.3829249254, 0.2417303085, 0.0605975520, 0.0062096803},
//     {0.0062096798, 0.0605975486, 0.2417303018, 0.3829249254, 0.2417303108, 0.0605975531, 0.0062096805},
//     {0.0062096796, 0.0605975475, 0.2417302996, 0.3829249254, 0.2417303130, 0.0605975542, 0.0062096806}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.0062096806, 0.0605975542, 0.2417303130, 0.3829249254, 0.2417302996, 0.0605975475, 0.0062096796},
//     {0.0062096805, 0.0605975531, 0.2417303108, 0.3829249254, 0.2417303018, 0.0605975486, 0.0062096798},
//     {0.0062096803, 0.0605975520, 0.2417303085, 0.3829249254, 0.2417303041, 0.0605975497, 0.0062096799},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096799, 0.0605975497, 0.2417303041, 0.3829249254, 0.2417303085, 0.0605975520, 0.0062096803},
//     {0.0062096798, 0.0605975486, 0.2417303018, 0.3829249254, 0.2417303108, 0.0605975531, 0.0062096805},
//     {0.0062096796, 0.0605975475, 0.2417302996, 0.3829249254, 0.2417303130, 0.0605975542, 0.0062096806}};

const double p_e = 0.000000001;
const double std_e = 0.01;

const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
const double risk_trans[7][7] = {
    {0.0062096802, 0.0605975512, 0.2417303070, 0.3829249254, 0.2417303056, 0.0605975505, 0.0062096801},
    {0.0062096801, 0.0605975511, 0.2417303067, 0.3829249254, 0.2417303059, 0.0605975506, 0.0062096801},
    {0.0062096801, 0.0605975510, 0.2417303065, 0.3829249254, 0.2417303061, 0.0605975508, 0.0062096801},
    {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
    {0.0062096801, 0.0605975508, 0.2417303061, 0.3829249254, 0.2417303065, 0.0605975510, 0.0062096801},
    {0.0062096801, 0.0605975506, 0.2417303059, 0.3829249254, 0.2417303067, 0.0605975511, 0.0062096801},
    {0.0062096801, 0.0605975505, 0.2417303056, 0.3829249254, 0.2417303070, 0.0605975512, 0.0062096802}};

const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
const double riskpre_trans[7][7] = {
    {0.0062096802, 0.0605975512, 0.2417303070, 0.3829249254, 0.2417303056, 0.0605975505, 0.0062096801},
    {0.0062096801, 0.0605975511, 0.2417303067, 0.3829249254, 0.2417303059, 0.0605975506, 0.0062096801},
    {0.0062096801, 0.0605975510, 0.2417303065, 0.3829249254, 0.2417303061, 0.0605975508, 0.0062096801},
    {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
    {0.0062096801, 0.0605975508, 0.2417303061, 0.3829249254, 0.2417303065, 0.0605975510, 0.0062096801},
    {0.0062096801, 0.0605975506, 0.2417303059, 0.3829249254, 0.2417303067, 0.0605975511, 0.0062096801},
    {0.0062096801, 0.0605975505, 0.2417303056, 0.3829249254, 0.2417303070, 0.0605975512, 0.0062096802}};

// const double risk_trans_preperiod[7][7] = {
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801}};

const double pi = 0.00005;

const double r_f = 0.03;

// Function Definitions:

// Utility
#define MUc(x) (pow((x), -rhopar))
#define inv_MU(u) (pow((u), (-(1.0 / rhopar))))
#define U(x) (pow((x), (1.0 - rhopar)) / (1.0 - rhopar))

// Grid
#define inter1d(x1, y1, y2) ((1.0 - (x1)) * (y1) + (x1) * (y2))
#define getwage(rrate) (1.0 - alphapar) * pow((alphapar / (rrate + deltapar)), (alphapar / (1.0 - alphapar)));
#define getlevel(x) (scale1 * (exp(exponen * (x)) + grmin))
#define getomega(x) ((x + 0.0) / size_j)
#define getgrid(x) (log((x) / scale1 - grmin) / exponen)

// EGM derivatives
#define nderiv(val1, val2, val3, x1, x2, x3) ((1.0 - (x3 - x2) / (x3 - x1)) * ((val3 - val2) / (x3 - x2)) + ((x3 - x2) / (x3 - x1)) * ((val2 - val1) / (x2 - x1)))

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

double K[size_k], Omega[size_j + 1];

void CreateFolder(const char *path)
{
    if (!CreateDirectory(path, NULL))
    {
        return;
    }
}

void copy(double *VectorIN, double *VectorOUT, int dim)
{
    int i;
    for (i = 0; i < dim; i++)
    {
        VectorOUT[i] = VectorIN[i];
    }
}

void null(double *VectorIN, int dim)
{
    int i;
    for (i = 0; i < dim; i++)
    {
        VectorIN[i] = 0;
    }
}

void POLICY(double *VF_final, double *dVF_final, double *save_final, double *VF, double *dVF, double *save, double *Portfolio, double K[size_k], double Omega[size_j], double wagerate)
{

    // INITIALIZATION //
    double *Kendo, *VFnew, *Kendo_min, temp, tempnext, dtempnext, *eVF, *deVF, critV, vfweight, slope1, slope2, tempvf, *consendo, *VFendo, *cohendo, cohexo, *VF_final_old;
    VFendo = (double *)calloc((ifulldim), sizeof(double)); // Value function on the next time grid, next iteration
    VFnew = (double *)calloc((ifulldim), sizeof(double));  // Value function on the next time grid, next iteration

    VF_final_old = (double *)calloc((iydim), sizeof(double)); // Value function on the next time grid, next iteration
    // Kendo = (double *)calloc((ifulldim), sizeof(double));  // endogenous grid values
    // Kendo_min = (double *)calloc((maxygrid), sizeof(double)); // endogenous grid values
    // eVF = (double *)calloc((ifulldim), sizeof(double));       // expected value function
    // deVF = (double *)calloc((ifulldim), sizeof(double));      // derivative of the expected value function
    cohendo = (double *)calloc((ifulldim), sizeof(double));  // Value function on the next time grid, next iteration
    consendo = (double *)calloc((ifulldim), sizeof(double)); // Value function on the next time grid, next iteration

    int i, ii, y, ynext, ypre, yprenext, j, iter, threshold_ii, Icase, itest, igridL, igridH, itemp;

    iter = 0;

    critV = 10000.0;
    // std::cout << "iter\t"
    //           << "critV\n";

    while (critV > epsV && iter < 300)
    {
        // we need copy to make a separate object

        null(cohendo, ifulldim);
        null(VFendo, ifulldim);

        copy(VF_final, VF_final_old, iydim);

        // std::cout << std::setprecision(16) << VF[inx(5, 5)] << "\n";
        // std::cout << std::setprecision(16) << VFnew[inx(5, 5)] << "\n";

        // main EGM computation
        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_y; y++)
            {
                for (ypre = 0; ypre < size_y; ypre++)
                {
                    // try omega here: omega index=j

                    for (j = 0; j < size_j + 1; j++)
                    {

                        tempnext = 0;
                        dtempnext = 0;

                        for (ynext = 0; ynext < size_y; ynext++)
                        {
                            tempnext += risk_trans[y][ynext] * VF[inx(i, ynext, y, j)];
                            dtempnext += risk_trans[y][ynext] * dVF[inx(i, ynext, y, j)];
                        }

                        cohendo[inx(i, y, ypre, j)] = K[i] + inv_MU(betapar * dtempnext);
                        VFendo[inx(i, y, ypre, j)] = U(cohendo[inx(i, y, ypre, j)] - K[i]) + betapar * tempnext;
                    }
                }
            }
        }

        // rescaling
        for (j = 0; j < size_j + 1; j++)
        {

            for (y = 0; y < size_y; y++)
            {
                for (ypre = 0; ypre < size_y; ypre++)
                {

                    threshold_ii = 0;

                    for (i = 0; i < size_k; i++)
                    {
                        // method 1: cash on hand
                        cohexo = (1.0 + (r_f + pi + risk_states[y]) * Omega[j] + r_f * (1 - Omega[j])) * K[i] + wagerate;

                        if (cohexo < cohendo[inx(0, y, ypre, j)])
                        {
                            save[inx(i, y, ypre, j)] = K[0];
                            VF[inx(i, y, ypre, j)] = U(cohexo - save[inx(i, y, ypre, j)]) + (VFendo[inx(0, y, ypre, j)] - U((cohendo[inx(0, y, ypre, j)] - K[0])));
                        }

                        if (cohexo >= cohendo[inx(0, y, ypre, j)])
                        {
                            itest = threshold_ii;

                            while ((itest < size_k) && cohexo > cohendo[(inx(itest, y, ypre, j))])
                            {
                                itest++;
                            }

                            if (itest == size_k)
                            {
                                // extrapolation
                                vfweight = (cohexo - cohendo[inx(size_k - 2, y, ypre, j)]) / (cohendo[inx(size_k - 1, y, ypre, j)] - cohendo[inx(size_k - 2, y, ypre, j)]);
                                igridL = size_k - 2;
                                igridH = size_k - 1;
                            }
                            else
                            {
                                // standard interior
                                vfweight = (cohexo - cohendo[inx(itest - 1, y, ypre, j)]) / (cohendo[inx(itest, y, ypre, j)] - cohendo[inx(itest - 1, y, ypre, j)]);
                                igridL = itest - 1;
                                igridH = itest - 0;
                            }

                            VF[inx(i, y, ypre, j)] = inter1d(vfweight, VFendo[inx(igridL, y, ypre, j)], VFendo[inx(igridH, y, ypre, j)]);
                            save[inx(i, y, ypre, j)] = inter1d(vfweight, K[igridL], K[igridH]);

                            threshold_ii = min(size_k - 2, itest);
                        }
                    }
                }
            }
        }

        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_y; y++)
            {
                for (ypre = 0; ypre < size_y; ypre++)
                {

                    for (j = 0; j < size_j + 1; j++)
                    {
                        tempnext = 0.0;
                        for (yprenext = 0; yprenext < size_y; yprenext++)
                        {
                            tempnext += riskpre_trans[ypre][yprenext] * VF[inx(i, yprenext, ypre, j)];
                        }

                        // std::cout << tempnext << "\n";
                        if (j == 0)
                        {
                            temp = tempnext;
                            itemp = 0;
                        }

                        if (tempnext > temp)
                        {
                            temp = tempnext;
                            itemp = j;
                        }
                    }

                    VF_final[inx2(i, y, ypre)] = VF[inx(i, y, ypre, itemp)];
                    save_final[inx2(i, y, ypre)] = save[inx(i, y, ypre, itemp)];
                    Portfolio[inx2(i, y, ypre)] = Omega[itemp];
                }
            }
        }

        // std::cout << std::setprecision(16) << VF[inx(5, 5)] << "\n";

        // computing new derivatives and convergence
        critV = 0.0;

        for (y = 0; y < size_y; y++)
        {
            for (ypre = 0; ypre < size_y; ypre++)
            {
                for (i = 0; i < size_k; i++)
                {

                    if (i >= 2)
                    {
                        dVF_final[inx2(i - 1, y, ypre)] = nderiv(VF_final[inx2(i - 2, y, ypre)], VF_final[inx2(i - 1, y, ypre)], VF_final[inx2(i, y, ypre)], K[i - 2], K[i - 1], K[i]);
                    }

                    critV = max(critV, abs(VF_final[inx2(i, y, ypre)] - VF_final_old[inx2(i, y, ypre)]));

                    // left corner
                    dVF_final[inx2(0, y, ypre)] = (VF_final[inx2(1, y, ypre)] - VF_final[inx2(0, y, ypre)]) / (K[1] - K[0]);
                    // right corner
                    dVF_final[inx2(size_k - 1, y, ypre)] = (VF_final[inx2(size_k - 1, y, ypre)] - VF_final[inx2(size_k - 2, y, ypre)]) / (K[size_k - 1] - K[size_k - 2]);
                }
            }
        }

        for (j = 0; j < size_j + 1; j++)
        {
            for (y = 0; y < size_y; y++)
            {
                for (ypre = 0; ypre < size_y; ypre++)
                {
                    for (i = 0; i < size_k; i++)
                    {
                        VF[inx(i, y, ypre, j)] = relaxVF * VF_final[inx2(i, y, ypre)] + (1 - relaxVF) * VF[inx(i, y, ypre, j)];
                    }
                }
            }
        }

        for (j = 0; j < size_j + 1; j++)
        {
            for (y = 0; y < size_y; y++)
            {
                for (ypre = 0; ypre < size_y; ypre++)
                {
                    for (i = 0; i < size_k; i++)
                    {

                        if (i >= 2)
                        {
                            dVF[inx(i - 1, y, ypre, j)] = nderiv(VF[inx(i - 2, y, ypre, j)], VF[inx(i - 1, y, ypre, j)], VF[inx(i, y, ypre, j)], K[i - 2], K[i - 1], K[i]);
                        }

                        // left corner
                        dVF[inx(0, y, ypre, j)] = (VF[inx(1, y, ypre, j)] - VF[inx(0, y, ypre, j)]) / (K[1] - K[0]);
                        // right corner
                        dVF[inx(size_k - 1, y, ypre, j)] = (VF[inx(size_k - 1, y, ypre, j)] - VF[inx(size_k - 2, y, ypre, j)]) / (K[size_k - 1] - K[size_k - 2]);
                    }
                }
            }
        }

        iter++;

        std::cout << "iteration=" << iter << ", critV=" << critV << "\n";
    }
}

void SIMULATION(double *save, double *dist, double *capitalout, double K[size_k])
{
    double *distold, critdist, distverif, weight;

    distold = (double *)calloc((iydim), sizeof(double));
    null(distold, iydim);

    int isave, i, y, ynext, ypre, yprenext, iter;

    critdist = 1.0;
    iter = 0;
    save[inx2(0, 0, 0)] = 0.001;
    while (critdist > epsdist && iter < 10000)
    {
        copy(dist, distold, iydim);
        null(dist, iydim);

        // std::cout << critdist << "\n";

        // distribution dynamics
        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_y; y++)
            {
                for (ypre = 0; ypre < size_y; ypre++)
                {
                    if (distold[inx2(i, y, ypre)] > 0)
                    {

                        isave = min((int)(floor(getgrid(save[inx2(i, y, ypre)]))), size_k - 2);
                        weight = (save[inx2(i, y, ypre)] - K[isave]) / (K[isave + 1] - K[isave]);
                        // std::cout << "weight=" << weight << "\n";

                        for (ynext = 0; ynext < size_y; ynext++)
                        {
                            for (yprenext = 0; yprenext < size_y; yprenext++)
                            {
                                dist[inx2(isave, ynext, yprenext)] += (1.0 - weight) * risk_trans[y][ynext] * riskpre_trans[ypre][yprenext] * distold[inx2(i, y, ypre)];
                                dist[inx2(min(isave + 1, size_k - 1), ynext, yprenext)] += (weight)*risk_trans[y][ynext] * riskpre_trans[ypre][yprenext] * (distold[inx2(i, y, ypre)]);
                            }
                        }
                    }
                }
            }
        }

        // convergence
        critdist = 0.0;
        distverif = 0.0;

        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_y; y++)
            {
                for (ypre = 0; ypre < size_y; ypre++)
                {
                    critdist = (max(critdist, abs(dist[inx2(i, y, ypre)] - distold[inx2(i, y, ypre)])));
                    distverif += dist[inx2(i, y, ypre)];
                }
            }
        }

        iter++;
        std::cout << "iteration time=" << iter << ", critdist=" << critdist << "\n";
        std::cout << "distverify=" << distverif << "\n";
    }

    *capitalout = 0.0;

    for (i = 0; i < size_k; i++)
    {
        for (y = 0; y < size_y; y++)
        {
            for (ypre = 0; ypre < size_y; ypre++)
            {
                *capitalout += dist[inx2(i, y, ypre)] * K[i];
            }
        }
    }
}

int main()
{
    // MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
    double *VF, *dVF, *save, *cons, *Portfolio;                                // for decision rules
    double capital1, capital0, PIB, critprice, taxL, welfare, rrate, wagerate; // for equilibrium

    // Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
    VF = (double *)calloc((ifulldim), sizeof(double));   // value function
    dVF = (double *)calloc((ifulldim), sizeof(double));  // value function derivative
    save = (double *)calloc((ifulldim), sizeof(double)); // value function

    double *VF_final, *dVF_final, *save_final;
    double *distin_final, *distout_final, capitalout; // for simulation

    VF_final = (double *)calloc((iydim), sizeof(double));  // value function
    dVF_final = (double *)calloc((iydim), sizeof(double)); // value function derivative

    save_final = (double *)calloc((iydim), sizeof(double));
    Portfolio = (double *)calloc((iydim), sizeof(double));
    // cons = (double *)calloc((ifulldim), sizeof(double));
    distin_final = (double *)calloc((iydim), sizeof(double));
    distout_final = (double *)calloc((iydim), sizeof(double));

    null(VF, ifulldim);
    null(dVF, ifulldim);

    null(VF_final, iydim);
    null(dVF_final, iydim);
    null(save_final, iydim);
    null(distin_final, iydim);
    null(distout_final, iydim);

    int i, y, ypre, j, tempcount;

    for (i = 0; i < size_k; i++)
    {
        K[i] = getlevel(i);
    }

    for (j = 0; j < size_j + 1; j++)
    {
        Omega[j] = getomega(j);
        // std::cout << Omega[j] << "\n";
    }

    // rrate = 0.040237086402090;
    wagerate = 0.6;
    distin_final[0] = 1.0;
    // taxL=0.3

    // initializing value function and initial derivatives

    for (j = 0; j < size_j + 1; j++)
    {
        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_y; y++)
            {
                for (ypre = 0; ypre < size_y; ypre++)
                {
                    VF[inx(i, y, ypre, j)] = U(wagerate + (1 + r_f + pi + risk_states[y]) * K[i]); // REQUIERE TO BE INCREASING IN K (the case here)
                    if (i >= 2)
                    {
                        dVF[inx(i - 1, y, ypre, j)] = nderiv(VF[inx(i - 2, y, ypre, j)], VF[inx(i - 1, y, ypre, j)], VF[inx(i, y, ypre, j)], K[i - 2], K[i - 1], K[i]);
                    }

                    // left corner
                    dVF[inx(0, y, ypre, j)] = (VF[inx(1, y, ypre, j)] - VF[inx(0, y, ypre, j)]) / (K[1] - K[0]);
                    // right corner
                    dVF[inx(size_k - 1, y, ypre, j)] = (VF[inx(size_k - 1, y, ypre, j)] - VF[inx(size_k - 2, y, ypre, j)]) / (K[size_k - 1] - K[size_k - 2]);
                }
            }
        }
    }

    for (i = 0; i < size_k; i++)
    {
        for (y = 0; y < size_y; y++)
        {
            for (ypre = 0; ypre < size_y; ypre++)
            {
                VF_final[inx2(i, y, ypre)] = U(wagerate + (1 + r_f + pi + risk_states[y]) * K[i]); // REQUIERE TO BE INCREASING IN K (the case here)
                if (i >= 2)
                {
                    dVF_final[inx2(i - 1, y, ypre)] = nderiv(VF_final[inx2(i - 2, y, ypre)], VF_final[inx2(i - 1, y, ypre)], VF_final[inx2(i, y, ypre)], K[i - 2], K[i - 1], K[i]);
                }

                // left corner
                dVF_final[inx2(0, y, ypre)] = (VF_final[inx2(1, y, ypre)] - VF_final[inx2(0, y, ypre)]) / (K[1] - K[0]);
                // right corner
                dVF_final[inx2(size_k - 1, y, ypre)] = (VF_final[inx2(size_k - 1, y, ypre)] - VF_final[inx2(size_k - 2, y, ypre)]) / (K[size_k - 1] - K[size_k - 2]);
            }
        }
    }

    POLICY(VF_final, dVF_final, save_final, VF, dVF, save, Portfolio, K, Omega, wagerate);
    printf("Policy Computation Done");
    // SIMULATION(save_final, distin_final, &capital1, K);

    // for (i = 0; i < size_k; i++)
    // {
    //     std::cout << i << "," << getlevel(i) << ",";

    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
    //     for (y = 0; y < size_y; y++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
    //         std::cout << distout[inx(i, y)] << ",";
    //     }
    //     std::cout << "\n";
    // }

    CreateFolder(".\\csv\\");
    CreateFolder(".\\figure\\");

    std::string filename_dist = "csv\\dist11,pe=e-9,std=0.01,premium=" + std::to_string(pi) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_j) + ",rho_c=" + std::to_string(rhopar) + ",Ksize=" + std::to_string(size_k) + ",relaxVF=" + std::to_string(relaxVF) + ".csv ";
    std::string filename_policy = "csv\\policy11,pe=e-9,std=0.01,premium=" + std::to_string(pi) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_j) + ",Ksize=" + std::to_string(size_k) + ",relaxVF=" + std::to_string(relaxVF) + ".csv ";
    std::string filename_VF = "csv\\VF11,pe=e-9,std=0.01,premium=" + std::to_string(pi) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_j) + ",Ksize=" + std::to_string(size_k) + ",relaxVF=" + std::to_string(relaxVF) + ".csv ";
    std::string filename_Port = "csv\\Portfolio11,pe=e-9,std=0.01,premium=" + std::to_string(pi) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_j) + ",Ksize=" + std::to_string(size_k) + ",relaxVF=" + std::to_string(relaxVF) + ".csv ";

    // std::string var = "sometext" + std::to_string(pi);
    // std::cout << var;

    // std::ofstream dfilecsv;
    // dfilecsv.open("csv\\dist4.csv");
    // dfilecsv << "gridnumber,"
    //          << "capital,"
    //          << "dist[0],"
    //          << "dist[1],"
    //          << "dist[2],"
    //          << "dist[3],"
    //          << "dist[4],"
    //          << "dist[5],"
    //          << "dist[6]\n";
    // for (i = 0; i < size_k; i++)
    // {
    //     dfilecsv << i << "," << getlevel(i) << ",";

    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
    //     for (y = 0; y < size_y; y++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
    //         if (y < size_y - 1)
    //         {
    //             dfilecsv << distin_final[inx2(i, y)] << ",";
    //         }

    //         if (y == size_y - 1)
    //         {
    //             dfilecsv << distin_final[inx2(i, y)];
    //         }
    //     }

    //     dfilecsv << "\n";
    // }

    // dfilecsv.close();

    std::ofstream policyfilecsv;
    policyfilecsv.open(filename_policy);

    policyfilecsv << "gridnumber,"
                  << "capital,";
    tempcount = 0;
    for (y = 0; y < size_y; y++)
    {
        for (ypre = 0; ypre < size_y; ypre++)
        {
            policyfilecsv << "policy[" << y << " " << ypre << "]";
            if (tempcount < size_y * size_y - 1)
            {
                policyfilecsv << ",";
            }
            else
            {
                policyfilecsv << "\n";
            }
            tempcount++;
        }
    }

    for (i = 0; i < size_k; i++)
    {
        policyfilecsv << i << "," << getlevel(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));

        for (y = 0; y < size_y; y++)

        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(i, ypre)]);
            for (ypre = 0; ypre < size_y; ypre++)
            {

                if (ypre + y * size_y < size_y * size_y - 1)
                {
                    policyfilecsv << save_final[inx2(i, y, ypre)] << ",";
                }

                if (ypre + y * size_y == size_y * size_y - 1)
                {
                    policyfilecsv << save_final[inx2(i, y, ypre)];
                }
            }
        }

        policyfilecsv << "\n";
    }

    policyfilecsv.close();

    // std::string var_VF = "csv\\VF,premium=" + std::to_string(premium) + ".csv";

    std::ofstream VFfilecsv;
    VFfilecsv.open(filename_VF);

    VFfilecsv << "gridnumber,"
              << "capital,";

    tempcount = 0;
    for (y = 0; y < size_y; y++)

    { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(i, ypre)]);
        for (ypre = 0; ypre < size_y; ypre++)
        {
            VFfilecsv << "VF[" << y << " " << ypre << "]";
            if (tempcount < size_y * size_y - 1)
            {
                VFfilecsv << ",";
            }
            else
            {
                VFfilecsv << "\n";
            }
            tempcount++;
        }
    }
    for (i = 0; i < size_k; i++)
    {
        VFfilecsv << i << "," << getlevel(i) << ",";
        for (y = 0; y < size_y; y++)

        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(i, ypre)]);
            for (ypre = 0; ypre < size_y; ypre++)
            {

                if (ypre + y * size_y < size_y * size_y - 1)
                {
                    VFfilecsv << VF_final[inx2(i, y, ypre)] << ",";
                }

                if (ypre + y * size_y == size_y * size_y - 1)
                {
                    VFfilecsv << VF_final[inx2(i, y, ypre)];
                }
            }
        }

        VFfilecsv << "\n";
    }

    VFfilecsv.close();

    // std::string var_port = "csv\\Portfolio,premium=" + std::to_string(premium) + ".csv";

    std::ofstream Portfilecsv;
    Portfilecsv.open(filename_Port);
    Portfilecsv << "gridnumber,"
                << "capital,";

    tempcount = 0;
    for (y = 0; y < size_y; y++)

    { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(i, ypre)]);
        for (ypre = 0; ypre < size_y; ypre++)
        {
            Portfilecsv << "Portfolio[" << y << " " << ypre << "]";

            if (tempcount < size_y * size_y - 1)
            {
                Portfilecsv << ",";
            }
            else
            {
                Portfilecsv << "\n";
            }
            tempcount++;
        }
    }

    for (i = 0; i < size_k; i++)
    {
        Portfilecsv << i << "," << getlevel(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < size_y; y++)

        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(i, ypre)]);
            for (ypre = 0; ypre < size_y; ypre++)
            {

                if (ypre + y * size_y < size_y * size_y - 1)
                {
                    Portfilecsv << Portfolio[inx2(i, y, ypre)] << ",";
                }

                if (ypre + y * size_y == size_y * size_y - 1)
                {
                    Portfilecsv << Portfolio[inx2(i, y, ypre)];
                }
            }
        }

        Portfilecsv << "\n";
    }

    Portfilecsv.close();
}
