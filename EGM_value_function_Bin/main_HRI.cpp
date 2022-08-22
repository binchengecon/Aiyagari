#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <direct.h>
#include <iomanip>
#include <cmath>
#include <limits>

const int size_k = 500; // number of grid points
const int size_m = 7;   // number of productivity classes
const int size_n = 7;
#define ifulldim (size_k * size_m * size_n)
#define inx(igridindex, jclassindex, kclassindex) (((kclassindex) * (size_m * size_k)) + ((jclassindex) * (size_k)) + (igridindex))

const double kmin = 0.0;
const double kmax = 500.0;

const double betapar = 0.9;
const double alphapar = 0.36;
const double deltapar = 0.08;
const double rhopar = 3.0;
const double labor = 1.0219882;

const double epsV = 1.0e-8;
const double epsdist = 1.0e-10;
const double epsK = 1.0e-6;
const double relaxsK = 0.01;

// grid constants
const double scale1 = 1.6;
const double grmin = (kmin / scale1) - 1.0;
const double exponen = log((kmax / scale1) - grmin) / (size_k - 1);

const double msstates[7] = {exp(-0.600000000000000), exp(-0.400000000000000), exp(-0.200000000000000), exp(0.000000000000000), exp(0.200000000000000), exp(0.400000000000000), exp(0.600000000000000)};

const double mstrans[7][7] = {
    {0.046746218637144, 0.217937777267117, 0.397822606398702, 0.266386738072197, 0.065169922261456, 0.005754191945237, 0.000182545418147},
    {0.023199661746751, 0.149524091076020, 0.369020347246402, 0.333823905199677, 0.110578117872631, 0.013276146769082, 0.000577730089437},
    {0.010548958644399, 0.093657511915497, 0.312761268311836, 0.382193227897354, 0.171253064028981, 0.027919224002876, 0.001666745199056},
    {0.004387354018187, 0.053538402796357, 0.242163972572887, 0.399820541225137, 0.242163972572887, 0.053538402796357, 0.004387354018187},
    {0.001666745199056, 0.027919224002876, 0.171253064028981, 0.382193227897354, 0.312761268311837, 0.093657511915497, 0.010548958644399},
    {0.000577730089436, 0.013276146769082, 0.110578117872631, 0.333823905199677, 0.369020347246403, 0.149524091076020, 0.023199661746751},
    {0.000182545418147, 0.005754191945237, 0.065169922261456, 0.266386738072197, 0.397822606398702, 0.217937777267117, 0.046746218637144}};

const double nsstates[7] = {exp(-0.600000000000000), exp(-0.400000000000000), exp(-0.200000000000000), exp(0.000000000000000), exp(0.200000000000000), exp(0.400000000000000), exp(0.600000000000000)};

const double nstrans[7][7] = {
    {0.046746218637144, 0.217937777267117, 0.397822606398702, 0.266386738072197, 0.065169922261456, 0.005754191945237, 0.000182545418147},
    {0.023199661746751, 0.149524091076020, 0.369020347246402, 0.333823905199677, 0.110578117872631, 0.013276146769082, 0.000577730089437},
    {0.010548958644399, 0.093657511915497, 0.312761268311836, 0.382193227897354, 0.171253064028981, 0.027919224002876, 0.001666745199056},
    {0.004387354018187, 0.053538402796357, 0.242163972572887, 0.399820541225137, 0.242163972572887, 0.053538402796357, 0.004387354018187},
    {0.001666745199056, 0.027919224002876, 0.171253064028981, 0.382193227897354, 0.312761268311837, 0.093657511915497, 0.010548958644399},
    {0.000577730089436, 0.013276146769082, 0.110578117872631, 0.333823905199677, 0.369020347246403, 0.149524091076020, 0.023199661746751},
    {0.000182545418147, 0.005754191945237, 0.065169922261456, 0.266386738072197, 0.397822606398702, 0.217937777267117, 0.046746218637144}};

// const double sstates[7] = {exp(-1.500000), exp(-1.000000), exp(-0.500000), exp(0.000000), exp(0.500000), exp(1.000000), exp(1.500000)};

// double strans[7][7] = {
//     {0.190787, 0.455383, 0.301749, 0.050061, 0.002002, 0.000019, 0.000000},
//     {0.052081, 0.301749, 0.455383, 0.173994, 0.016424, 0.000367, 0.000002},
//     {0.008774, 0.121520, 0.419444, 0.365696, 0.080233, 0.004279, 0.000053},
//     {0.000889, 0.029507, 0.235589, 0.468029, 0.235589, 0.029507, 0.000889},
//     {0.000053, 0.004279, 0.080233, 0.365696, 0.419444, 0.121520, 0.008774},
//     {0.000002, 0.000367, 0.016424, 0.173994, 0.455383, 0.301749, 0.052081},
//     {0.000000, 0.000019, 0.002002, 0.050061, 0.301749, 0.455383, 0.190787}};

// Functions:

// Utility
#define MUc(x) (pow((x), -rhopar))
#define inv_MU(u) (pow((u), (-(1.0 / rhopar))))
#define U(x) (pow((x), (1.0 - rhopar)) / (1.0 - rhopar))

// Grid
#define inter1d(x1, y1, y2) ((1.0 - (x1)) * (y1) + (x1) * (y2))
#define getwage(rrate) (1.0 - alphapar) * pow((alphapar / (rrate + deltapar)), (alphapar / (1.0 - alphapar)));
#define getlevel(x) (scale1 * (exp(exponen * (x)) + grmin))
#define getgrid(x) (log((x) / scale1 - grmin) / exponen)

// EGM derivatives
#define nderiv(val1, val2, val3, x1, x2, x3) ((1.0 - (x3 - x2) / (x3 - x1)) * ((val3 - val2) / (x3 - x2)) + ((x3 - x2) / (x3 - x1)) * ((val2 - val1) / (x2 - x1)))

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

double K[size_k];

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

void POLICY_WF(double *VF, double *dVF, double *save, double wagerate, double rrate, double K[size_k])
{

    // INITIALIZATION //
    double *Kendo, *VFnew, *Kendo_min, tempnext, dtempnext, *eVF, *deVF, critV, vfweight, slope1, slope2, tempvf, *consendo, *VFendo, *cohendo, cohexo;
    VFendo = (double *)calloc((ifulldim), sizeof(double)); // Value function on the next time grid, next iteration
    VFnew = (double *)calloc((ifulldim), sizeof(double));  // Value function on the next time grid, next iteration
    // Kendo = (double *)calloc((ifulldim), sizeof(double));  // endogenous grid values
    // Kendo_min = (double *)calloc((maxygrid), sizeof(double)); // endogenous grid values
    // eVF = (double *)calloc((ifulldim), sizeof(double));       // expected value function
    // deVF = (double *)calloc((ifulldim), sizeof(double));      // derivative of the expected value function
    cohendo = (double *)calloc((ifulldim), sizeof(double));  // cash on hand on the next time grid, next iteration
    consendo = (double *)calloc((ifulldim), sizeof(double)); // consumption on the next time grid, next iteration

    int i, ii, y, ynext, n, nnext, iter, threshold_ii, Icase, itest, igridL, igridH;

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

        // std::cout << std::setprecision(16) << VF[inx(5, 5)] << "\n";
        // std::cout << std::setprecision(16) << VFnew[inx(5, 5)] << "\n";

        // main EGM computation
        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_m; y++)
            {
                for (n = 0; n < size_n; n++)
                {

                    tempnext = 0;
                    dtempnext = 0;

                    for (ynext = 0; ynext < size_m; ynext++)
                    {
                        for (nnext = 0; nnext < size_n; nnext++)
                        {
                            tempnext += mstrans[y][ynext] * nstrans[n][nnext] * VF[inx(i, ynext, nnext)];
                            dtempnext += mstrans[y][ynext] * nstrans[n][nnext] * dVF[inx(i, ynext, nnext)];
                        }
                    }

                    cohendo[inx(i, y, n)] = K[i] + inv_MU(betapar * dtempnext);
                    VFendo[inx(i, y, n)] = U(cohendo[inx(i, y, n)] - K[i]) + betapar * tempnext;
                }
            }
        }

        // rescaling
        for (y = 0; y < size_m; y++)
        {
            for (n = 0; n < size_n; n++)
            {
                threshold_ii = 0;

                for (i = 0; i < size_k; i++)
                {
                    // method 1: cash on hand
                    cohexo = (1.0 + nsstates[n] * rrate) * K[i] + wagerate * msstates[y];

                    if (cohexo < cohendo[inx(0, y, n)])
                    {
                        save[inx(i, y, n)] = K[0];
                        VF[inx(i, y, n)] = U(cohexo - save[inx(i, y, n)]) + (VFendo[inx(0, y, n)] - U((cohendo[inx(0, y, n)] - K[0])));
                    }

                    if (cohexo >= cohendo[inx(0, y, n)])
                    {
                        itest = threshold_ii;

                        while ((itest < size_k) && cohexo > cohendo[(inx(itest, y, n))])
                        {
                            itest++;
                        }

                        // if (itest>size_k){
                        //     // error
                        // }

                        if (itest == size_k)
                        {
                            // extrapolation
                            vfweight = (cohexo - cohendo[inx(size_k - 2, y, n)]) / (cohendo[inx(size_k - 1, y, n)] - cohendo[inx(size_k - 2, y, n)]);
                            igridL = size_k - 2;
                            igridH = size_k - 1;
                        }
                        else
                        {
                            // standard interior
                            vfweight = (cohexo - cohendo[inx(itest - 1, y, n)]) / (cohendo[inx(itest, y, n)] - cohendo[inx(itest - 1, y, n)]);
                            igridL = itest - 1;
                            igridH = itest - 0;
                        }

                        VF[inx(i, y, n)] = inter1d(vfweight, VFendo[inx(igridL, y, n)], VFendo[inx(igridH, y, n)]);
                        save[inx(i, y, n)] = inter1d(vfweight, K[igridL], K[igridH]);

                        threshold_ii = min(size_k - 2, itest);
                    }
                }
            }
        }

        // std::cout << std::setprecision(16) << VF[inx(5, 5)] << "\n";

        // computing new derivatives and convergence
        critV = 0.0;

        for (y = 0; y < size_m; y++)
        {
            for (n = 0; n < size_n; n++)
            {
                for (i = 0; i < size_k; i++)
                {

                    if (i >= 2)
                    {
                        dVF[inx(i - 1, y, n)] = nderiv(VF[inx(i - 2, y, n)], VF[inx(i - 1, y, n)], VF[inx(i, y, n)], K[i - 2], K[i - 1], K[i]);
                    }

                    critV = max(critV, abs(VF[inx(i, y, n)] - VFnew[inx(i, y, n)]));

                    // left corner
                    dVF[inx(0, y, n)] = (VF[inx(1, y, n)] - VF[inx(0, y, n)]) / (K[1] - K[0]);
                    // right corner
                    dVF[inx(size_k - 1, y, n)] = (VF[inx(size_k - 1, y, n)] - VF[inx(size_k - 2, y, n)]) / (K[size_k - 1] - K[size_k - 2]);
                }
            }
        }

        // std::cout << "check derivative numbers: start\n";
        // std::cout << std::setprecision(16) << VF[inx(4, 5)] << "\n";
        // std::cout << std::setprecision(16) << VF[inx(5, 5)] << "\n";
        // std::cout << std::setprecision(16) << VF[inx(6, 5)] << "\n";

        // std::cout << std::setprecision(16) << K[4] << "\n";
        // std::cout << std::setprecision(16) << K[5] << "\n";
        // std::cout << std::setprecision(16) << K[6] << "\n";
        // std::cout << "check derivative numbers: end\n";

        // std::cout << std::setprecision(16) << dVF[inx(5, 5)] << "\n";

        // std::cout << std::setprecision(16) << dVF[inx(5, 5)] << "\n";

        iter++;

        // std::cout << iter << "\t" << std::setprecision(15) << critV << "\n";
    }
}

void SIMULATION_WF(double *save, double *dist, double *capitalout, double K[size_k])
{
    double *distold, critdist, distverif, weight;
    distold = (double *)calloc((ifulldim), sizeof(double));
    null(distold, ifulldim);

    int isave, i, y, ynext, n, nnext;

    critdist = 1.0;

    while (critdist > epsdist)
    {
        copy(dist, distold, ifulldim);
        null(dist, ifulldim);

        // distribution dynamics
        for (i = 0; i < size_k; i++)
        {
            for (y = 0; y < size_m; y++)
            {

                for (n = 0; n < size_n; n++)
                {
                    if (distold[inx(i, y, n)] > 0)
                    {

                        isave = min((int)(floor(getgrid(save[inx(i, y, n)]))), size_k - 2);
                        weight = (save[inx(i, y, n)] - K[isave]) / (K[isave + 1] - K[isave]);

                        for (ynext = 0; ynext < size_m; ynext++)
                        {
                            for (nnext = 0; nnext < size_n; nnext++)
                            {
                                dist[inx(isave, ynext, nnext)] += (1.0 - weight) * mstrans[y][ynext] * nstrans[n][nnext] * distold[inx(i, y, nnext)];
                                dist[inx(min(isave + 1, size_k - 1), ynext, nnext)] += (weight)*mstrans[y][ynext] * nstrans[n][nnext] * (distold[inx(i, y, nnext)]);
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
            for (y = 0; y < size_m; y++)
            {
                for (n = 0; n < size_n; n++)
                {
                    critdist = (max(critdist, abs(dist[inx(i, y, n)] - distold[inx(i, y, n)])));
                    distverif += dist[inx(i, y, n)];
                }
            }
        }

        // std::cout << critdist << "\n";
    }

    *capitalout = 0.0;

    for (i = 0; i < size_k; i++)
    {
        for (y = 0; y < size_m; y++)
        {
            for (n = 0; n < size_n; n++)
            {
                *capitalout += dist[inx(i, y, n)] * K[i];
            }
        }
    }
}

int main()
{
    // MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
    double *VF, *dVF, *save, *cons;                                            // for decision rules
    double *distin, *distout, capitalout;                                      // for simulation
    double capital1, capital0, PIB, critprice, taxL, welfare, rrate, wagerate; // for equilibrium

    // Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
    VF = (double *)calloc((ifulldim), sizeof(double));  // value function
    dVF = (double *)calloc((ifulldim), sizeof(double)); // value function
    save = (double *)calloc((ifulldim), sizeof(double));
    // cons = (double *)calloc((ifulldim), sizeof(double));
    distin = (double *)calloc((ifulldim), sizeof(double));
    distout = (double *)calloc((ifulldim), sizeof(double));

    null(VF, ifulldim);
    null(dVF, ifulldim);
    null(save, ifulldim);
    null(distin, ifulldim);
    null(distout, ifulldim);

    int i, y, n;

    for (i = 0; i < size_k; i++)
    {
        K[i] = getlevel(i);
    }

    rrate = 0.040237086402090;
    wagerate = getwage(rrate);
    distin[0] = 1.0;
    // taxL=0.3

    // initializing value function and initial derivatives

    for (i = 0; i < size_k; i++)
    {
        for (y = 0; y < size_m; y++)
        {
            for (n = 0; n < size_n; n++)
            {
                VF[inx(i, y, n)] = U(msstates[y] * wagerate + (1 + nsstates[n] * rrate) * K[i]); // REQUIERE TO BE INCREASING IN K (the case here)

                if (i >= 2)
                {
                    dVF[inx(i - 1, y, n)] = nderiv(VF[inx(i - 2, y, n)], VF[inx(i - 1, y, n)], VF[inx(i, y, n)], K[i - 2], K[i - 1], K[i]);
                }

                // left corner
                dVF[inx(0, y, n)] = (VF[inx(1, y, n)] - VF[inx(0, y, n)]) / (K[1] - K[0]);
                // right corner
                dVF[inx(size_k - 1, y, n)] = (VF[inx(size_k - 1, y, n)] - VF[inx(size_k - 2, y, n)]) / (K[size_k - 1] - K[size_k - 2]);
            }
        }
    }

    // // No Firms
    // POLICY(VF, dVF, save, cons, wagerate, rrate, K);
    // std::cout << "policy computation is over\n";
    // SIMULATION(save, distin, distout, K);
    // std::cout << "simulation computation is over\n";

    // With Firms
    double critK = 1.0;
    double capitalin, rrate0;

    std::cout << "critK\t"
              << "r_in\t"
              << "r_out\t"
              << "K_in\t"
              << "K_out\t"
              << "L\t"
              << "\n";

    while (critK > epsK)
    {

        wagerate = getwage(rrate);

        capitalin = labor * pow((alphapar / (rrate + deltapar)), (1.0 / (1.0 - alphapar)));

        POLICY_WF(VF, dVF, save, wagerate, rrate, K);

        SIMULATION_WF(save, distin, &capital1, K);

        rrate0 = rrate;

        rrate = (alphapar * (pow(labor, (1.0 - alphapar)))) / (pow(((relaxsK * capital1 + (1.0 - relaxsK) * capitalin)), (1.0 - alphapar))) - deltapar;

        critK = abs((capital1 - capitalin) / capitalin);

        std::cout << critK << "\t" << rrate0 << "\t" << rrate << "\t" << capitalin << "\t" << capital1 << "\t" << labor << "\t"
                  << "\n";
    }

    // for (i = 0; i < ifulldim; i++)
    // {
    //     std::cout
    //         << save[i] - K[i] << "\t";
    // }

    // for (i = 0; i < size_k; i++)
    // {
    //     std::cout << i << "," << getlevel(i) << ",";

    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
    //     for (y = 0; y < size_m; y++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
    //         std::cout << distout[inx(i, y)] << ",";
    //     }
    //     std::cout << "\n";
    // }

    std::ofstream dfilecsv;
    dfilecsv.open("csv\\dist_HRI.csv");
    dfilecsv << "gridnumber,"
             << "capital,\n";

    for (i = 0; i < size_k; i++)
    {
        dfilecsv << i << "," << getlevel(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < size_m; y++)
        { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
            for (n = 0; n < size_n; n++)
            {
                dfilecsv << distin[inx(i, y, n)] << ",";
            }
        }
        dfilecsv << "\n";
    }

    dfilecsv.close();

    // std::ofstream policyfilecsv;
    // policyfilecsv.open("csv\\policy_HRI.csv");
    // policyfilecsv << "gridnumber,"
    //               << "capital,"
    //               << "policy[0],"
    //               << "policy[1],"
    //               << "policy[2],"
    //               << "policy[3],"
    //               << "policy[4],"
    //               << "policy[5],"
    //               << "policy[6],\n";
    // for (i = 0; i < size_k; i++)
    // {
    //     policyfilecsv << i << "," << getlevel(i) << ",";

    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
    //     for (y = 0; y < size_m; y++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
    //         policyfilecsv << save[inx(i, y)] << ",";
    //     }
    //     policyfilecsv << "\n";
    // }

    // policyfilecsv.close();

    std::ofstream VFfilecsv;
    VFfilecsv.open("csv\\VF_HRI.csv");
    VFfilecsv << "gridnumber,"
              << "capital,\n";
    for (i = 0; i < size_k; i++)
    {
        VFfilecsv << i << "," << getlevel(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < size_m; y++)
        { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
            for (n = 0; n < size_n; n++)
            {
                VFfilecsv << VF[inx(i, y, n)] << ",";
            }
        }
        VFfilecsv << "\n";
    }

    VFfilecsv.close();
}
