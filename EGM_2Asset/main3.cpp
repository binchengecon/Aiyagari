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

const int size_k = 500; // number of grid points
const int size_y = 7;   // number of productivity classes
const int size_j = 100;

#define ifulldim (size_k * size_y * size_j)
#define iydim (size_k * size_y)
#define inx(igridindex, jclassindex, kclassindex) (((kclassindex) * (size_k * size_y)) + ((jclassindex) * (size_k)) + (igridindex))
#define inx2(igridindex, jclassindex) (((jclassindex) * (size_k)) + (igridindex))

const double kmin = 0.0;
const double kmax = 500.0;

const double betapar = 0.96;
const double alphapar = 0.36;
const double deltapar = 0.08;
const double rhopar = 3.0;
const double labor = 1.0219882;

const double epsV = 1.0e-8;
const double epsdist = 1.0e-10;
const double epsK = 1.0e-6;
const double relaxsK = 0.005;

// grid constants
const double scale1 = 1.6;
const double grmin = (kmin / scale1) - 1.0;
const double exponen = log((kmax / scale1) - grmin) / (size_k - 1);

const double sstates[7] = {-0.030619, -0.020412, -0.010206, 0.000000, 0.010206, 0.020412, 0.030619};
const double strans[7][7] = {
    {0.026240, 0.152924, 0.361483, 0.328567, 0.114742, 0.015266, 0.000778},
    {0.016044, 0.114742, 0.328567, 0.361483, 0.152924, 0.024700, 0.001539},
    {0.009452, 0.082835, 0.287445, 0.382789, 0.196114, 0.038437, 0.002929},
    {0.005362, 0.057531, 0.242024, 0.390166, 0.242024, 0.057531, 0.005362},
    {0.002929, 0.038437, 0.196114, 0.382789, 0.287445, 0.082835, 0.009452},
    {0.001539, 0.024700, 0.152924, 0.361483, 0.328567, 0.114742, 0.016044},
    {0.000778, 0.015266, 0.114742, 0.328567, 0.361483, 0.152924, 0.026240}};

const double pi = 0.05;

const double r_f = 0.02;

// Function Definitions:

// Utility
#define MUc(x) (pow((x), -rhopar))
#define inv_MU(u) (pow((u), (-(1.0 / rhopar))))
#define U(x) (pow((x), (1.0 - rhopar)) / (1.0 - rhopar))

// Grid
#define inter1d(x1, y1, y2) ((1.0 - (x1)) * (y1) + (x1) * (y2))
#define getwage(rrate) (1.0 - alphapar) * pow((alphapar / (rrate + deltapar)), (alphapar / (1.0 - alphapar)));
#define getlevel(x) (scale1 * (exp(exponen * (x)) + grmin))
#define getomega(x) ((x + 1.0) / size_j)
#define getgrid(x) (log((x) / scale1 - grmin) / exponen)

// EGM derivatives
#define nderiv(val1, val2, val3, x1, x2, x3) ((1.0 - (x3 - x2) / (x3 - x1)) * ((val3 - val2) / (x3 - x2)) + ((x3 - x2) / (x3 - x1)) * ((val2 - val1) / (x2 - x1)))

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

double K[size_k], Omega[size_j];

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

void POLICY(double *VF_final, double *dVF_final, double *save_final, double *VF, double *dVF, double *save, double K[size_k], double Omega[size_j], double wagerate)
{

    // INITIALIZATION //
    double *Kendo, *VFnew, *Kendo_min, tempnext, dtempnext, *eVF, *deVF, critV, vfweight, slope1, slope2, tempvf, *consendo, *VFendo, *cohendo, cohexo, *VF_final_old;
    VFendo = (double *)calloc((ifulldim), sizeof(double)); // Value function on the next time grid, next iteration
    VFnew = (double *)calloc((ifulldim), sizeof(double));  // Value function on the next time grid, next iteration

    VF_final_old = (double *)calloc((iydim), sizeof(double)); // Value function on the next time grid, next iteration
    // Kendo = (double *)calloc((ifulldim), sizeof(double));  // endogenous grid values
    // Kendo_min = (double *)calloc((maxygrid), sizeof(double)); // endogenous grid values
    // eVF = (double *)calloc((ifulldim), sizeof(double));       // expected value function
    // deVF = (double *)calloc((ifulldim), sizeof(double));      // derivative of the expected value function
    cohendo = (double *)calloc((ifulldim), sizeof(double));  // Value function on the next time grid, next iteration
    consendo = (double *)calloc((ifulldim), sizeof(double)); // Value function on the next time grid, next iteration

    int i, ii, y, ynext, j, iter, threshold_ii, Icase, itest, igridL, igridH;

    iter = 0;

    critV = 10000.0;
    // std::cout << "iter\t"
    //           << "critV\n";

    while (critV > epsV)
    {
        // we need copy to make a separate object
        copy(VF, VFnew, ifulldim);

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

                // try omega here: omega index=j

                for (j = 0; j < size_j; j++)
                {

                    tempnext = 0;
                    dtempnext = 0;

                    for (ynext = 0; ynext < size_y; ynext++)
                    {
                        tempnext += strans[y][ynext] * VF[inx(i, ynext, j)];
                        dtempnext += strans[y][ynext] * dVF[inx(i, ynext, j)];
                    }

                    cohendo[inx(i, y, j)] = K[i] + inv_MU(betapar * dtempnext);
                    VFendo[inx(i, y, j)] = U(cohendo[inx(i, y, j)] - K[i]) + betapar * tempnext;
                }
            }
        }

        // rescaling
        for (j = 0; j < size_j; j++)
        {

            for (y = 0; y < size_y; y++)
            {

                threshold_ii = 0;

                for (i = 0; i < size_k; i++)
                {
                    // method 1: cash on hand
                    cohexo = (1.0 + (r_f + pi + sstates[y]) * Omega[j] + r_f * (1 - Omega[j])) * K[i] + wagerate;

                    if (cohexo < cohendo[inx(0, y, j)])
                    {
                        save[inx(i, y, j)] = K[0];
                        VF[inx(i, y, j)] = U(cohexo - save[inx(i, y, j)]) + (VFendo[inx(0, y, j)] - U((cohendo[inx(0, y, j)] - K[0])));
                    }

                    if (cohexo >= cohendo[inx(0, y, j)])
                    {
                        itest = threshold_ii;

                        while ((itest < size_k) && cohexo > cohendo[(inx(itest, y, j))])
                        {
                            itest++;
                        }

                        if (itest == size_k)
                        {
                            // extrapolation
                            vfweight = (cohexo - cohendo[inx(size_k - 2, y, j)]) / (cohendo[inx(size_k - 1, y, j)] - cohendo[inx(size_k - 2, y, j)]);
                            igridL = size_k - 2;
                            igridH = size_k - 1;
                        }
                        else
                        {
                            // standard interior
                            vfweight = (cohexo - cohendo[inx(itest - 1, y, j)]) / (cohendo[inx(itest, y, j)] - cohendo[inx(itest - 1, y, j)]);
                            igridL = itest - 1;
                            igridH = itest - 0;
                        }

                        VF[inx(i, y, j)] = inter1d(vfweight, VFendo[inx(igridL, y, j)], VFendo[inx(igridH, y, j)]);
                        save[inx(i, y, j)] = inter1d(vfweight, K[igridL], K[igridH]);

                        threshold_ii = min(size_k - 2, itest);
                    }
                }
            }
        }

        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_y; y++)
            {

                double temp;
                int itemp;

                for (j = 0; j < size_j; j++)
                {
                    tempnext = 0.0;
                    for (ynext = 0; ynext < size_y; ynext++)
                    {
                        tempnext += strans[y][ynext] * VF[inx(i, ynext, j)];
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

                VF_final[inx2(i, y)] = VF[inx(i, y, itemp)];
                save_final[inx2(i, y)] = save[inx(i, y, itemp)];
            }
        }

        // std::cout << std::setprecision(16) << VF[inx(5, 5)] << "\n";

        // computing new derivatives and convergence
        critV = 0.0;

        for (y = 0; y < size_y; y++)
        {
            for (i = 0; i < size_k; i++)
            {

                if (i >= 2)
                {
                    dVF_final[inx2(i - 1, y)] = nderiv(VF_final[inx2(i - 2, y)], VF_final[inx2(i - 1, y)], VF_final[inx2(i, y)], K[i - 2], K[i - 1], K[i]);
                }

                critV = max(critV, abs(VF_final[inx2(i, y)] - VF_final_old[inx2(i, y)]));

                // left corner
                dVF_final[inx2(0, y)] = (VF_final[inx2(1, y)] - VF_final[inx2(0, y)]) / (K[1] - K[0]);
                // right corner
                dVF_final[inx2(size_k - 1, y)] = (VF_final[inx2(size_k - 1, y)] - VF_final[inx2(size_k - 2, y)]) / (K[size_k - 1] - K[size_k - 2]);
            }
        }

        for (j = 0; j < size_j; j++)
        {
            for (y = 0; y < size_y; y++)
            {
                for (i = 0; i < size_k; i++)
                {

                    if (i >= 2)
                    {
                        dVF[inx(i - 1, y, j)] = nderiv(VF_final[inx2(i - 2, y)], VF_final[inx2(i - 1, y)], VF_final[inx2(i, y)], K[i - 2], K[i - 1], K[i]);
                    }

                    // left corner
                    dVF[inx(0, y, j)] = (VF_final[inx2(1, y)] - VF_final[inx2(0, y)]) / (K[1] - K[0]);
                    // right corner
                    dVF[inx(size_k - 1, y, j)] = (VF_final[inx2(size_k - 1, y)] - VF_final[inx2(size_k - 2, y)]) / (K[size_k - 1] - K[size_k - 2]);
                }
            }
        }

        iter++;

        std::cout << critV << "\n";
    }
}

void SIMULATION(double *save, double *dist, double *capitalout, double K[size_k])
{
    double *distold, critdist, distverif, weight;

    distold = (double *)calloc((iydim), sizeof(double));
    null(distold, iydim);

    int isave, i, y, ynext;

    critdist = 1.0;

    while (critdist > epsdist)
    {
        copy(dist, distold, iydim);
        null(dist, iydim);

        std::cout << critdist << "\n";

        // distribution dynamics
        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_y; y++)
            {

                if (distold[inx2(i, y)] > 0)
                {

                    isave = min((int)(floor(getgrid(save[inx2(i, y)]))), size_k - 2);
                    weight = (save[inx2(i, y)] - K[isave]) / (K[isave + 1] - K[isave]);

                    for (ynext = 0; ynext < size_y; ynext++)
                    {
                        dist[inx2(isave, ynext)] += (1.0 - weight) * strans[y][ynext] * distold[inx2(i, y)];
                        dist[inx2(min(isave + 1, size_k - 1), ynext)] += (weight)*strans[y][ynext] * (distold[inx2(i, y)]);
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
                critdist = (max(critdist, abs(dist[inx2(i, y)] - distold[inx2(i, y)])));
                distverif += dist[inx2(i, y)];
            }
        }

        // std::cout << critdist << "\n";
    }

    *capitalout = 0.0;

    for (i = 0; i < size_k; i++)
    {
        for (y = 0; y < size_y; y++)
        {
            *capitalout += dist[inx2(i, y)] * K[i];
        }
    }
}

int main()
{
    // MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
    double *VF, *dVF, *save, *cons;                                            // for decision rules
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

    int i, y, j;

    for (i = 0; i < size_k; i++)
    {
        K[i] = getlevel(i);
    }

    for (j = 0; j < size_j; j++)
    {
        Omega[j] = getomega(j);
        // std::cout << Omega[j] << "\n";
    }

    rrate = 0.040237086402090;
    wagerate = 0.03;
    distin_final[0] = 1.0;
    // taxL=0.3

    // initializing value function and initial derivatives

    for (j = 0; j < size_j; j++)
    {
        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_y; y++)
            {
                VF[inx(i, y, j)] = U(wagerate + (1 + rrate) * K[i]); // REQUIERE TO BE INCREASING IN K (the case here)
                if (i >= 2)
                {
                    dVF[inx(i - 1, y, j)] = nderiv(VF[inx(i - 2, y, j)], VF[inx(i - 1, y, j)], VF[inx(i, y, j)], K[i - 2], K[i - 1], K[i]);
                }

                // left corner
                dVF[inx(0, y, j)] = (VF[inx(1, y, j)] - VF[inx(0, y, j)]) / (K[1] - K[0]);
                // right corner
                dVF[inx(size_k - 1, y, j)] = (VF[inx(size_k - 1, y, j)] - VF[inx(size_k - 2, y, j)]) / (K[size_k - 1] - K[size_k - 2]);
            }
        }
    }

    for (i = 0; i < size_k; i++)
    {
        for (y = 0; y < size_y; y++)
        {
            VF_final[inx2(i, y)] = U(wagerate + (1 + rrate) * K[i]); // REQUIERE TO BE INCREASING IN K (the case here)
            if (i >= 2)
            {
                dVF_final[inx2(i - 1, y)] = nderiv(VF_final[inx2(i - 2, y)], VF_final[inx2(i - 1, y)], VF_final[inx2(i, y)], K[i - 2], K[i - 1], K[i]);
            }

            // left corner
            dVF_final[inx2(0, y)] = (VF_final[inx2(1, y)] - VF_final[inx2(0, y)]) / (K[1] - K[0]);
            // right corner
            dVF_final[inx2(size_k - 1, y)] = (VF_final[inx2(size_k - 1, y)] - VF_final[inx2(size_k - 2, y)]) / (K[size_k - 1] - K[size_k - 2]);
        }
    }

    POLICY(VF_final, dVF_final, save_final, VF, dVF, save, K, Omega, wagerate);

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

    // std::ofstream dfilecsv;
    // dfilecsv.open("csv\\dist.csv");
    // dfilecsv << "gridnumber,"
    //          << "capital,"
    //          << "dist[0],"
    //          << "dist[1],"
    //          << "dist[2],"
    //          << "dist[3],"
    //          << "dist[4],"
    //          << "dist[5],"
    //          << "dist[6],\n";
    // for (i = 0; i < size_k; i++)
    // {
    //     dfilecsv << i << "," << getlevel(i) << ",";

    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
    //     for (y = 0; y < size_y; y++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
    //         dfilecsv << distin_final[inx2(i, y)] << ",";
    //     }
    //     dfilecsv << "\n";
    // }

    // dfilecsv.close();

    CreateFolder(".\\csv\\");
    CreateFolder(".\\figure\\");

    std::ofstream policyfilecsv;
    policyfilecsv.open("csv\\policy.csv");
    policyfilecsv << "gridnumber,"
                  << "capital,"
                  << "policy[0],"
                  << "policy[1],"
                  << "policy[2],"
                  << "policy[3],"
                  << "policy[4],"
                  << "policy[5],"
                  << "policy[6]\n";
    for (i = 0; i < size_k; i++)
    {
        policyfilecsv << i << "," << getlevel(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < size_y; y++)
        { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
            if (y < size_y - 1)
            {
                policyfilecsv << save_final[inx2(i, y)] << ",";
            }

            if (y == size_y - 1)
            {
                policyfilecsv << save_final[inx2(i, y)];
            }
        }
        policyfilecsv << "\n";
    }

    policyfilecsv.close();

    std::ofstream VFfilecsv;
    VFfilecsv.open("csv\\VF.csv");
    VFfilecsv << "gridnumber,"
              << "capital,"
              << "VF[0],"
              << "VF[1],"
              << "VF[2],"
              << "VF[3],"
              << "VF[4],"
              << "VF[5],"
              << "VF[6]\n";
    for (i = 0; i < size_k; i++)
    {
        VFfilecsv << i << "," << getlevel(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < size_y; y++)
        { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
            if (y < size_y - 1)
            {
                VFfilecsv << VF_final[inx2(i, y)] << ",";
            }

            if (y == size_y - 1)
            {
                VFfilecsv << VF_final[inx2(i, y)];
            }
        }
        VFfilecsv << "\n";
    }

    VFfilecsv.close();
}
