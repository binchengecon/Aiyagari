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

const int size_asset = 500; // number of grid points
const int size_risk = 7;    // number of productivity classes
const int size_laborincome = 7;
const int size_portfoliochoice = 100;

#define ARLP_dim (size_asset * size_risk * size_laborincome * (size_portfoliochoice + 1))
#define ARL_dim (size_asset * size_risk * size_laborincome)
#define index_ARLP(asset_gridindex, risk_gridindex, laborincome_gridindex, portfoliochoice_gridindex) (((portfoliochoice_gridindex) * (size_asset * size_risk * size_laborincome)) + ((laborincome_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))

#define index_ARL(asset_gridindex, risk_gridindex, laborincome_gridindex) (((laborincome_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))

const double kmin = 0.0;
const double kmax = 500.0;

const double betapar = 0.9;
const double alphapar = 0.36;
const double deltapar = 0.08;
const double rhopar = 3.0;
const double labor = 1.0219882;

const double epsV = 1.0e-8;
const double epsdist = 1.0e-9;
const double epsK = 1.0e-6;
const double relaxsK = 0.005;

// grid constants
const double scale1 = 1.6;
const double grmin = (kmin / scale1) - 1.0;
const double exponen = log((kmax / scale1) - grmin) / (size_asset - 1);

const double laborincome_states[7] = {exp(-0.600000000000000), exp(-0.400000000000000), exp(-0.200000000000000), exp(0.000000000000000), exp(0.200000000000000), exp(0.400000000000000), exp(0.600000000000000)};

const double laborincome_trans[7][7] = {
    {0.046746218637144, 0.217937777267117, 0.397822606398702, 0.266386738072197, 0.065169922261456, 0.005754191945237, 0.000182545418147},
    {0.023199661746751, 0.149524091076020, 0.369020347246402, 0.333823905199677, 0.110578117872631, 0.013276146769082, 0.000577730089437},
    {0.010548958644399, 0.093657511915497, 0.312761268311836, 0.382193227897354, 0.171253064028981, 0.027919224002876, 0.001666745199056},
    {0.004387354018187, 0.053538402796357, 0.242163972572887, 0.399820541225137, 0.242163972572887, 0.053538402796357, 0.004387354018187},
    {0.001666745199056, 0.027919224002876, 0.171253064028981, 0.382193227897354, 0.312761268311837, 0.093657511915497, 0.010548958644399},
    {0.000577730089436, 0.013276146769082, 0.110578117872631, 0.333823905199677, 0.369020347246403, 0.149524091076020, 0.023199661746751},
    {0.000182545418147, 0.005754191945237, 0.065169922261456, 0.266386738072197, 0.397822606398702, 0.217937777267117, 0.046746218637144}};

// p_e = 0.2
// const double risk_states[7] = {-0.030619, -0.020412, -0.010206, 0.000000, 0.010206, 0.020412, 0.030619};
// const double risk_trans[7][7] = {
//     {0.026240, 0.152924, 0.361483, 0.328567, 0.114742, 0.015266, 0.000778},
//     {0.016044, 0.114742, 0.328567, 0.361483, 0.152924, 0.024700, 0.001539},
//     {0.009452, 0.082835, 0.287445, 0.382789, 0.196114, 0.038437, 0.002929},
//     {0.005362, 0.057531, 0.242024, 0.390166, 0.242024, 0.057531, 0.005362},
//     {0.002929, 0.038437, 0.196114, 0.382789, 0.287445, 0.082835, 0.009452},
//     {0.001539, 0.024700, 0.152924, 0.361483, 0.328567, 0.114742, 0.016044},
//     {0.000778, 0.015266, 0.114742, 0.328567, 0.361483, 0.152924, 0.026240}};

//  p_e=0.001
const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
const double risk_trans[7][7] = {
    {0.006262, 0.060934, 0.242398, 0.382924, 0.241063, 0.060262, 0.006157},
    {0.006245, 0.060822, 0.242175, 0.382924, 0.241285, 0.060374, 0.006175},
    {0.006227, 0.060710, 0.241953, 0.382925, 0.241508, 0.060486, 0.006192},
    {0.006210, 0.060597, 0.241730, 0.382925, 0.241730, 0.060597, 0.006210},
    {0.006192, 0.060486, 0.241508, 0.382925, 0.241953, 0.060710, 0.006227},
    {0.006175, 0.060374, 0.241285, 0.382924, 0.242175, 0.060822, 0.006245},
    {0.006157, 0.060262, 0.241063, 0.382924, 0.242398, 0.060934, 0.006262}};

const double premium = 0.001; // 0.0005, 0.0001 not converging, stuck at 3.34268e-006????
// const double premium = 0.0000;

const double r_f = 0.03;
// const double r_f = 0.040237086402090;

// Function Definitions:

// Utility
#define MUc(x) (pow((x), -rhopar))
#define inv_MU(u) (pow((u), (-(1.0 / rhopar))))
#define U(x) (pow((x), (1.0 - rhopar)) / (1.0 - rhopar))

// Grid
#define inter1d(x1, y1, y2) ((1.0 - (x1)) * (y1) + (x1) * (y2))
#define getwage(rrate) (1.0 - alphapar) * pow((alphapar / (rrate + deltapar)), (alphapar / (1.0 - alphapar)));
#define getlevel(x) (scale1 * (exp(exponen * (x)) + grmin))
#define getomega(x) ((x + 0.0) / size_portfoliochoice)
#define getgrid(x) (log((x) / scale1 - grmin) / exponen)

// EGM derivatives
#define nderiv(val1, val2, val3, x1, x2, x3) ((1.0 - (x3 - x2) / (x3 - x1)) * ((val3 - val2) / (x3 - x2)) + ((x3 - x2) / (x3 - x1)) * ((val2 - val1) / (x2 - x1)))

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

double K[size_asset],
    Omega[size_portfoliochoice + 1];

void CreateFolder(const char *path)
{
    if (!CreateDirectory(path, NULL))
    {
        return;
    }
}

void copy(double *VectorIN, double *VectorOUT, int dim)
{
    int asset_index;
    for (asset_index = 0; asset_index < dim; asset_index++)
    {
        VectorOUT[asset_index] = VectorIN[asset_index];
    }
}

void null(double *VectorIN, int dim)
{
    int asset_index;
    for (asset_index = 0; asset_index < dim; asset_index++)
    {
        VectorIN[asset_index] = 0;
    }
}

void POLICY(double *VF, double *dVF, double *save, double *Portfolio, double K[size_asset], double Omega[size_portfoliochoice + 1], double wagerate)
{

    // INITIALIZATION //
    double *Kendo, *VFnew, *Kendo_min, *OmegaEndo, tempnext, dtempnext, *eVF, *deVF, critV, vfweight, weight, slope1, slope2, temp, dtemp, *consendo, *VFendo, *cohendo, cohexo, *VF_old, next_period_wealth;
    VFendo = (double *)calloc((ARL_dim), sizeof(double)); // Value function on the next time grid, next iteration
    OmegaEndo = (double *)calloc((ARL_dim), sizeof(double));
    VF_old = (double *)calloc((ARL_dim), sizeof(double)); // Value function on the next time grid, next iteration
    // Kendo = (double *)calloc((ARLP_dim), sizeof(double));  // endogenous grid values
    // Kendo_min = (double *)calloc((maxygrid), sizeof(double)); // endogenous grid values
    // eVF = (double *)calloc((ARLP_dim), sizeof(double));       // expected value function
    // deVF = (double *)calloc((ARLP_dim), sizeof(double));      // derivative of the expected value function
    cohendo = (double *)calloc((ARL_dim), sizeof(double));  // Value function on the next time grid, next iteration
    consendo = (double *)calloc((ARL_dim), sizeof(double)); // Value function on the next time grid, next iteration

    int asset_index, ii, risk_index, risk_indexnext, laborincome_index, laborincome_indexnext, portfoliochoice_index, iter, threshold_ii, Icase, itest, igridL, igridH, isavenext, itemp;
    // std::cout<< "PlicyComputationStart";
    iter = 0;

    critV = 10000.0;
    // std::cout << "iter\t"
    //           << "critV\n";

    while (critV > epsV)
    {
        // we need copy to make a separate object
        // copy(VF, VFnew, ARL_dim);

        null(cohendo, ARL_dim);
        null(VFendo, ARL_dim);

        copy(VF, VF_old, ARL_dim);

        // std::cout << std::setprecision(16) << VF[index_ARLP(5, 5)] << "\n";
        // std::cout << std::setprecision(16) << VFnew[index_ARLP(5, 5)] << "\n";

        // main EGM computation
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {

                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {

                    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
                    {

                        tempnext = 0;
                        dtempnext = 0;

                        next_period_wealth = K[asset_index] * (1.0 + r_f + (premium + risk_states[risk_index]) * Omega[portfoliochoice_index]);

                        // here find
                        isavenext = min((int)(floor(getgrid(next_period_wealth))), size_asset - 1);
                        weight = (next_period_wealth - K[isavenext]) / (K[isavenext + 1] - K[isavenext]);

                        for (risk_indexnext = 0; risk_indexnext < size_risk; risk_indexnext++)
                        {
                            for (laborincome_indexnext = 0; laborincome_indexnext < size_laborincome; laborincome_indexnext++)
                            {
                                tempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * inter1d(weight, VF[index_ARL(isavenext, risk_indexnext, laborincome_indexnext)], VF[index_ARL(isavenext + 1, risk_indexnext, laborincome_indexnext)]);
                                dtempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * inter1d(weight, dVF[index_ARL(isavenext, risk_indexnext, laborincome_indexnext)], dVF[index_ARL(isavenext + 1, risk_indexnext, laborincome_indexnext)]);
                            }
                        }

                        // here need to find the best tempnext
                        if (portfoliochoice_index == 0)
                        {
                            temp = tempnext;
                            dtemp = dtempnext;
                            itemp = 0;
                        }

                        if (tempnext > temp)
                        {
                            temp = tempnext;
                            dtemp = dtempnext;
                            itemp = portfoliochoice_index;
                        }
                    }

                    cohendo[index_ARL(asset_index, risk_index, laborincome_index)] = K[asset_index] + inv_MU(betapar * dtemp);
                    VFendo[index_ARL(asset_index, risk_index, laborincome_index)] = U(cohendo[index_ARL(asset_index, risk_index, laborincome_index)] - K[asset_index]) + betapar * temp;
                    OmegaEndo[index_ARL(asset_index, risk_index, laborincome_index)] = Omega[itemp];
                }
            }

            // try omega here: omega index=portfoliochoice_index
        }

        // rescaling
        // for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
        // {

        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                threshold_ii = 0;

                for (asset_index = 0; asset_index < size_asset; asset_index++)
                {
                    // method 1: cash on hand
                    cohexo = K[asset_index] + wagerate * laborincome_states[laborincome_index];

                    if (cohexo < cohendo[index_ARL(0, risk_index, laborincome_index)])
                    {
                        save[index_ARL(asset_index, risk_index, laborincome_index)] = K[0];
                        VF[index_ARL(asset_index, risk_index, laborincome_index)] = U(cohexo - save[index_ARL(asset_index, risk_index, laborincome_index)]) + (VFendo[index_ARL(0, risk_index, laborincome_index)] - U((cohendo[index_ARL(0, risk_index, laborincome_index)] - K[0])));
                        Portfolio[index_ARL(asset_index, risk_index, laborincome_index)] = (OmegaEndo[index_ARL(0, risk_index, laborincome_index)]);
                    }

                    if (cohexo >= cohendo[index_ARL(0, risk_index, laborincome_index)])
                    {
                        itest = threshold_ii;

                        while ((itest < size_asset) && cohexo > cohendo[(index_ARL(itest, risk_index, laborincome_index))])
                        {
                            itest++;
                        }

                        if (itest == size_asset)
                        {
                            // extrapolation
                            vfweight = (cohexo - cohendo[index_ARL(size_asset - 2, risk_index, laborincome_index)]) / (cohendo[index_ARL(size_asset - 1, risk_index, laborincome_index)] - cohendo[index_ARL(size_asset - 2, risk_index, laborincome_index)]);
                            igridL = size_asset - 2;
                            igridH = size_asset - 1;
                        }
                        else
                        {
                            // standard interior
                            vfweight = (cohexo - cohendo[index_ARL(itest - 1, risk_index, laborincome_index)]) / (cohendo[index_ARL(itest, risk_index, laborincome_index)] - cohendo[index_ARL(itest - 1, risk_index, laborincome_index)]);
                            igridL = itest - 1;
                            igridH = itest - 0;
                        }

                        VF[index_ARL(asset_index, risk_index, laborincome_index)] = inter1d(vfweight, VFendo[index_ARL(igridL, risk_index, laborincome_index)], VFendo[index_ARL(igridH, risk_index, laborincome_index)]);
                        save[index_ARL(asset_index, risk_index, laborincome_index)] = inter1d(vfweight, K[igridL], K[igridH]);

                        Portfolio[index_ARL(asset_index, risk_index, laborincome_index)] = inter1d(vfweight, OmegaEndo[index_ARL(igridL, risk_index, laborincome_index)], OmegaEndo[index_ARL(igridH, risk_index, laborincome_index)]);

                        threshold_ii = min(size_asset - 2, itest);
                    }
                }
            }
        } // }

        // std::cout << std::setprecision(16) << VF[index_ARLP(5, 5)] << "\n";

        // computing new derivatives and convergence
        critV = 0.0;

        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                for (asset_index = 0; asset_index < size_asset; asset_index++)
                {

                    if (asset_index >= 2)
                    {
                        dVF[index_ARL(asset_index - 1, risk_index, laborincome_index)] = nderiv(VF[index_ARL(asset_index - 2, risk_index, laborincome_index)], VF[index_ARL(asset_index - 1, risk_index, laborincome_index)], VF[index_ARL(asset_index, risk_index, laborincome_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                    }

                    critV = max(critV, abs(VF[index_ARL(asset_index, risk_index, laborincome_index)] - VF_old[index_ARL(asset_index, risk_index, laborincome_index)]));

                    // left corner
                    dVF[index_ARL(0, risk_index, laborincome_index)] = (VF[index_ARL(1, risk_index, laborincome_index)] - VF[index_ARL(0, risk_index, laborincome_index)]) / (K[1] - K[0]);
                    // right corner
                    dVF[index_ARL(size_asset - 1, risk_index, laborincome_index)] = (VF[index_ARL(size_asset - 1, risk_index, laborincome_index)] - VF[index_ARL(size_asset - 2, risk_index, laborincome_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                }
            }
        }

        iter++;

        std::cout << critV << "\n";
    }
}

void SIMULATION(double *save, double *dist, double *capitalout, double K[size_asset])
{
    double *distold, critdist, distverif, weight;

    distold = (double *)calloc((ARL_dim), sizeof(double));
    null(distold, ARL_dim);

    int isave, asset_index, risk_index, risk_indexnext, laborincome_index, laborincome_indexnext, iter;

    critdist = 1.0;
    iter = 0;

    while (critdist > epsdist && iter < 3000)
    {
        copy(dist, distold, ARL_dim);
        null(dist, ARL_dim);

        // std::cout << critdist << "\n";

        // distribution dynamics
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    if (distold[index_ARL(asset_index, risk_index, laborincome_index)] > 0)
                    {

                        isave = min((int)(floor(getgrid(save[index_ARL(asset_index, risk_index, laborincome_index)]))), size_asset - 1);
                        weight = (save[index_ARL(asset_index, risk_index, laborincome_index)] - K[isave]) / (K[isave + 1] - K[isave]);
                        // std::cout << "weight=" << weight << "\n";

                        for (risk_indexnext = 0; risk_indexnext < size_risk; risk_indexnext++)
                        {
                            for (laborincome_indexnext = 0; laborincome_indexnext < size_laborincome; laborincome_indexnext++)
                            {
                                dist[index_ARL(isave, risk_indexnext, laborincome_indexnext)] += (1.0 - weight) * risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * distold[index_ARL(asset_index, risk_index, laborincome_index)];
                                dist[index_ARL(min(isave + 1, size_asset - 1), risk_indexnext, laborincome_indexnext)] += (weight)*risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * (distold[index_ARL(asset_index, risk_index, laborincome_index)]);
                            }
                        }
                    }
                }
            }
        }

        // convergence
        critdist = 0.0;
        distverif = 0.0;

        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    critdist = (max(critdist, abs(dist[index_ARL(asset_index, risk_index, laborincome_index)] - distold[index_ARL(asset_index, risk_index, laborincome_index)])));
                    distverif += dist[index_ARL(asset_index, risk_index, laborincome_index)];
                }
            }
        }

        iter++;
    }
    std::cout << "iteration=" << iter << ", critdist=" << critdist << ", distverify=" << distverif << "\n";

    *capitalout = 0.0;

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                *capitalout += dist[index_ARL(asset_index, risk_index, laborincome_index)] * K[asset_index];
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
    VF = (double *)calloc((ARL_dim), sizeof(double));   // value function
    dVF = (double *)calloc((ARL_dim), sizeof(double));  // value function derivative
    save = (double *)calloc((ARL_dim), sizeof(double)); // value function

    double *VF_final, *dVF_final, *save_final;
    double *distin_final, *distout_final, capitalout; // for simulation

    Portfolio = (double *)calloc((ARL_dim), sizeof(double));
    // cons = (double *)calloc((ARLP_dim), sizeof(double));
    distin_final = (double *)calloc((ARL_dim), sizeof(double));
    distout_final = (double *)calloc((ARL_dim), sizeof(double));

    null(VF, ARL_dim);
    null(dVF, ARL_dim);

    null(save, ARL_dim);
    null(distin_final, ARL_dim);
    null(distout_final, ARL_dim);

    int asset_index, risk_index, laborincome_index, portfoliochoice_index, tempcount;

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        K[asset_index] = getlevel(asset_index);
    }

    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    {
        // Omega[portfoliochoice_index] = getomega(portfoliochoice_index);
        Omega[portfoliochoice_index] = 0;
        // std::cout << Omega[portfoliochoice_index] << "\n";
    }

    rrate = 0.040237086402090;
    // wagerate = 1.18;
    wagerate = 0.34;
    distin_final[0] = 1.0;
    // taxL=0.3

    // initializing value function and initial derivatives

    // for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    // {
    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                VF[index_ARL(asset_index, risk_index, laborincome_index)] = U(wagerate + (1 + rrate) * K[asset_index]); // REQUIERE TO BE INCREASING IN K (the case here)
            }
        }
    }
    // }

    // for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    // {
    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                if (asset_index >= 2)
                {
                    dVF[index_ARL(asset_index - 1, risk_index, laborincome_index)] = nderiv(VF[index_ARL(asset_index - 2, risk_index, laborincome_index)], VF[index_ARL(asset_index - 1, risk_index, laborincome_index)], VF[index_ARL(asset_index, risk_index, laborincome_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                }

                // left corner
                dVF[index_ARL(0, risk_index, laborincome_index)] = (VF[index_ARL(1, risk_index, laborincome_index)] - VF[index_ARL(0, risk_index, laborincome_index)]) / (K[1] - K[0]);
                // right corner
                dVF[index_ARL(size_asset - 1, risk_index, laborincome_index)] = (VF[index_ARL(size_asset - 1, risk_index, laborincome_index)] - VF[index_ARL(size_asset - 2, risk_index, laborincome_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
            }
        }
    }
    // }

    // for (asset_index = 0; asset_index < size_asset; asset_index++)
    // {
    //     for (risk_index = 0; risk_index < size_risk; risk_index++)
    //     {
    //         for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    //         {
    //             VF_final[index_ARL(asset_index, risk_index, laborincome_index)] = U(wagerate + (1 + rrate) * K[asset_index]); // REQUIERE TO BE INCREASING IN K (the case here)
    //         }
    //     }
    // }

    // for (asset_index = 0; asset_index < size_asset; asset_index++)
    // {
    //     for (risk_index = 0; risk_index < size_risk; risk_index++)
    //     {
    //         for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    //         {
    //             if (asset_index >= 2)
    //             {
    //                 dVF_final[index_ARL(asset_index - 1, risk_index, laborincome_index)] = nderiv(VF_final[index_ARL(asset_index - 2, risk_index, laborincome_index)], VF_final[index_ARL(asset_index - 1, risk_index, laborincome_index)], VF_final[index_ARL(asset_index, risk_index, laborincome_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
    //             }

    //             // left corner
    //             dVF_final[index_ARL(0, risk_index, laborincome_index)] = (VF_final[index_ARL(1, risk_index, laborincome_index)] - VF_final[index_ARL(0, risk_index, laborincome_index)]) / (K[1] - K[0]);
    //             // right corner
    //             dVF_final[index_ARL(size_asset - 1, risk_index, laborincome_index)] = (VF_final[index_ARL(size_asset - 1, risk_index, laborincome_index)] - VF_final[index_ARL(size_asset - 2, risk_index, laborincome_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
    //         }
    //     }
    // }

    POLICY(VF, dVF, save, Portfolio, K, Omega, wagerate);
    printf("Policy Computation Done\n");
    SIMULATION(save, distin_final, &capital1, K);

    // for (asset_index = 0; asset_index < size_asset; asset_index++)
    // {
    //     std::cout << asset_index << "," << getlevel(asset_index) << ",";

    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
    //     for (risk_index = 0; risk_index < size_risk; risk_index++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(asset_index, risk_index)]);
    //         std::cout << distout[index_ARLP(asset_index, risk_index)] << ",";
    //     }
    //     std::cout << "\n";
    // }

    CreateFolder(".\\csv\\");
    CreateFolder(".\\figure\\");

    std::string filename_dist = "csv\\dist8,premium=" + std::to_string(premium) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_portfoliochoice) + ".csv ";
    std::string filename_policy = "csv\\policy8,premium=" + std::to_string(premium) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_portfoliochoice) + ".csv ";
    std::string filename_VF = "csv\\VF8,premium=" + std::to_string(premium) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_portfoliochoice) + ".csv ";
    std::string filename_Port = "csv\\Portfolio8,premium=" + std::to_string(premium) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_portfoliochoice) + ".csv ";

    std::ofstream dfilecsv;
    dfilecsv.open(filename_dist);
    // dfilecsv << "gridnumber,"
    //          << "capital,"
    //          << "dist[0],"
    //          << "dist[1],"
    //          << "dist[2],"
    //          << "dist[3],"
    //          << "dist[4],"
    //          << "dist[5],"
    //          << "dist[6]\n";
    dfilecsv << "gridnumber,"
             << "capital,";
    tempcount = 0;
    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            dfilecsv << "dist[" << laborincome_index << " " << risk_index << "]";
            if (tempcount < size_laborincome * size_risk - 1)
            {
                dfilecsv << ",";
            }
            else
            {
                dfilecsv << "\n";
            }
            tempcount++;
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        dfilecsv << asset_index << "," << getlevel(asset_index) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(asset_index, risk_index)]);
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {

                if (risk_index + laborincome_index * size_risk < size_risk * size_laborincome - 1)
                {
                    dfilecsv << distin_final[index_ARL(asset_index, risk_index, laborincome_index)] << ",";
                }

                if (risk_index + laborincome_index * size_risk == size_risk * size_laborincome - 1)
                {
                    dfilecsv << distin_final[index_ARL(asset_index, risk_index, laborincome_index)];
                }
            }
        }

        dfilecsv << "\n";
    }

    dfilecsv.close();

    // std::string var_policy = "csv\\policy,premium=" + std::to_string(premium) + ".csv";

    std::ofstream policyfilecsv;
    policyfilecsv.open(filename_policy);

    policyfilecsv << "gridnumber,"
                  << "capital,";
    tempcount = 0;
    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            policyfilecsv << "policy[" << laborincome_index << " " << risk_index << "]";
            if (tempcount < size_laborincome * size_risk - 1)
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

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        policyfilecsv << asset_index << "," << getlevel(asset_index) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));

        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(asset_index, risk_index)]);
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {

                if (risk_index + laborincome_index * size_risk < size_risk * size_laborincome - 1)
                {
                    policyfilecsv << save[index_ARL(asset_index, risk_index, laborincome_index)] << ",";
                }

                if (risk_index + laborincome_index * size_risk == size_risk * size_laborincome - 1)
                {
                    policyfilecsv << save[index_ARL(asset_index, risk_index, laborincome_index)];
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
    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            VFfilecsv << "VF[" << laborincome_index << " " << risk_index << "]";
            if (tempcount < size_laborincome * size_risk - 1)
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
    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        VFfilecsv << asset_index << "," << getlevel(asset_index) << ",";
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(asset_index, risk_index)]);
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {

                if (risk_index + laborincome_index * size_risk < size_risk * size_laborincome - 1)
                {
                    VFfilecsv << VF[index_ARL(asset_index, risk_index, laborincome_index)] << ",";
                }

                if (risk_index + laborincome_index * size_risk == size_risk * size_laborincome - 1)
                {
                    VFfilecsv << VF[index_ARL(asset_index, risk_index, laborincome_index)];
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
    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
        {
            Portfilecsv << "Portfolio[" << laborincome_index << " " << risk_index << "]";

            if (tempcount < size_laborincome * size_risk - 1)
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

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        Portfilecsv << asset_index << "," << getlevel(asset_index) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(asset_index, risk_index)]);
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {

                if (risk_index + laborincome_index * size_risk < size_risk * size_laborincome - 1)
                {
                    Portfilecsv << Portfolio[index_ARL(asset_index, risk_index, laborincome_index)] << ",";
                }

                if (risk_index + laborincome_index * size_risk == size_risk * size_laborincome - 1)
                {
                    Portfilecsv << Portfolio[index_ARL(asset_index, risk_index, laborincome_index)];
                }
            }
        }

        Portfilecsv << "\n";
    }

    Portfilecsv.close();
}
