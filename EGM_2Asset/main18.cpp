// Modify the correlation of risky return and labor income

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <sys/time.h>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <iostream>
// #include <direct.h>
#include <iomanip>
#include <cmath>
#include <limits>
// #include <windows.h>
#include <string>
#include <sstream>
#include "txt/Shock.hpp"
#include "txt/Risk_Labor.hpp"

// const int size_risk = 1; // number of productivity classes
// const int size_shock = 7;
// const int size_laborincome = 7;

const int size_asset = 20; // number of grid points
const int size_portfoliochoice = 4;

#define ARRLLRRLLP_dim (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock * size_shock * (size_portfoliochoice + 1))
#define ARRLLRRLL_dim (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock * size_shock)

#define index_ARRLLRRLLP(asset_gridindex, risk_gridindex, risk_pre_gridindex, laborincome_gridindex, laborincome_pre_gridindex, riskshock_gridindex, riskshock_pre_gridindex, laborshock_gridindex, laborshock_pre_gridindex, portfoliochoice_gridindex) (((portfoliochoice_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock * size_shock)) + ((laborshock_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock)) + ((laborshock_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock)) + ((riskshock_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock)) + ((riskshock_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome)) + ((laborincome_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome)) + ((laborincome_gridindex) * (size_asset * size_risk * size_risk)) + ((risk_pre_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))
#define index_ARRLLRRLL(asset_gridindex, risk_gridindex, risk_pre_gridindex, laborincome_gridindex, laborincome_pre_gridindex, riskshock_gridindex, riskshock_pre_gridindex, laborshock_gridindex, laborshock_pre_gridindex) (((laborshock_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock)) + ((laborshock_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock)) + ((riskshock_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock)) + ((riskshock_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome)) + ((laborincome_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome)) + ((laborincome_gridindex) * (size_asset * size_risk * size_risk)) + ((risk_pre_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))

const double kmin = 0.0;
const double kmax = 300.0;

const double betapar = 0.9;
const double alphapar = 0.36;
const double deltapar = 0.08;
const double rhopar = 3.0;
const double labor = 1.0219882;

const double epsV = 1.0e-8;
const double epsdist = 1.0e-8;
const double epsK = 1.0e-6;
const double relaxsK = 0.005;
const double relaxVF = 0.000;

// grid constants
const double scale1 = 1.6;
const double grmin = (kmin / scale1) - 1.0;
const double exponen = log((kmax / scale1) - grmin) / (size_asset - 1);

const double pi = 0.00005;
const double corr = 0.8;
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
#define getomega(x) ((x + 0.0) / size_portfoliochoice)
#define getgrid(x) (log((x) / scale1 - grmin) / exponen)

// EGM derivatives
#define nderiv(val1, val2, val3, x1, x2, x3) ((1.0 - (x3 - x2) / (x3 - x1)) * ((val3 - val2) / (x3 - x2)) + ((x3 - x2) / (x3 - x1)) * ((val2 - val1) / (x2 - x1)))

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

double K[size_asset], Omega[size_portfoliochoice + 1];

// void CreateFolder(const char *path)
// {
//     if (!CreateDirectory(path, NULL))
//     {
//         return;
//     }
// }

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

void POLICY(double *VF_final, double *dVF_final, double *save_final, double *VF, double *dVF, double *save, double *Portfolio, double K[size_asset], double Omega[size_portfoliochoice + 1], double wagerate)
{

    // INITIALIZATION //
    double *Kendo, *VFnew, *Kendo_min, temp, tempnext, dtempnext, *eVF, *deVF, critV, vfweight, slope1, slope2, tempvf, *consendo, *VFendo, *cohendo, cohexo, *VF_final_old, conditional_measure;
    VFendo = (double *)calloc((ARRLLRRLLP_dim), sizeof(double)); // Value function on the next time grid, next iteration
    VFnew = (double *)calloc((ARRLLRRLLP_dim), sizeof(double));  // Value function on the next time grid, next iteration

    VF_final_old = (double *)calloc((ARRLLRRLL_dim), sizeof(double)); // Value function on the next time grid, next iteration
    // Kendo = (double *)calloc((ARRLLRRLLP_dim), sizeof(double));  // endogenous grid values
    // Kendo_min = (double *)calloc((maxygrid), sizeof(double)); // endogenous grid values
    // eVF = (double *)calloc((ARRLLRRLLP_dim), sizeof(double));       // expected value function
    // deVF = (double *)calloc((ARRLLRRLLP_dim), sizeof(double));      // derivative of the expected value function
    cohendo = (double *)calloc((ARRLLRRLLP_dim), sizeof(double));  // Value function on the next time grid, next iteration
    consendo = (double *)calloc((ARRLLRRLLP_dim), sizeof(double)); // Value function on the next time grid, next iteration

    int shockstate_pre, shockstate_current, shockstate_next, asset_index, ii, risk_index, risk_indexnext, risk_pre_index, risk_pre_indexnext, laborincome_index, laborincome_indexnext, laborincome_pre_index, laborincome_pre_indexnext, riskshock_index, riskshock_indexnext, riskshock_pre_index, riskshock_pre_indexnext, laborshock_index, laborshock_indexnext, laborshock_pre_index, laborshock_pre_indexnext, portfoliochoice_index, iter, threshold_ii, Icase, itest, igridL, igridH, itemp;

    iter = 0;

    critV = 10000.0;
    // std::cout << "iter\t"
    //           << "critV\n";

    while (critV > epsV && iter < 250)
    {
        // we need copy to make a separate object

        null(cohendo, ARRLLRRLLP_dim);
        null(VFendo, ARRLLRRLLP_dim);

        copy(VF_final, VF_final_old, ARRLLRRLL_dim);

        // std::cout << std::setprecision(16) << VF[index_ARRLLRRLLP(5, 5)] << "\n";
        // std::cout << std::setprecision(16) << VFnew[index_ARRLLRRLLP(5, 5)] << "\n";

        // main EGM computation
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    // try omega here: omega index=portfoliochoice_index
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {

                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
                                            {

                                                tempnext = 0;
                                                dtempnext = 0;
                                                shockstate_current = riskshock_index * size_shock + laborshock_index;
                                                for (risk_indexnext = 0; risk_indexnext < size_risk; risk_indexnext++)
                                                {
                                                    for (laborincome_indexnext = 0; laborincome_indexnext < size_laborincome; laborincome_indexnext++)
                                                    {
                                                        for (riskshock_indexnext = 0; riskshock_indexnext < size_shock; riskshock_indexnext++)
                                                        {
                                                            for (laborshock_indexnext = 0; laborshock_indexnext < size_shock; laborshock_indexnext++)
                                                            {
                                                                shockstate_next = riskshock_indexnext * size_shock + laborshock_indexnext;

                                                                tempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * risk_labor_shock_trans[0][shockstate_next] * VF[index_ARRLLRRLLP(asset_index, risk_indexnext, risk_index, laborincome_indexnext, laborincome_index, riskshock_indexnext, riskshock_index, laborshock_indexnext, laborshock_index, portfoliochoice_index)];
                                                                dtempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * risk_labor_shock_trans[0][shockstate_next] * dVF[index_ARRLLRRLLP(asset_index, risk_indexnext, risk_index, laborincome_indexnext, laborincome_index, riskshock_indexnext, riskshock_index, laborshock_indexnext, laborshock_index, portfoliochoice_index)];
                                                            }
                                                        }
                                                    }
                                                }

                                                cohendo[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = K[asset_index] + inv_MU(betapar * dtempnext);
                                                VFendo[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = U(cohendo[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - K[asset_index]) + betapar * tempnext;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // std::cout << "EGM done\n";

        // rescalings
        for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
        {

            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            shockstate_current = riskshock_index * size_shock + laborshock_index;

                                            threshold_ii = 0;

                                            for (asset_index = 0; asset_index < size_asset; asset_index++)
                                            {
                                                // method 1: cash on hand
                                                cohexo = (1.0 + (r_f + pi + risk_states[risk_index]) * risk_shock_states[shockstate_current] * Omega[portfoliochoice_index] + r_f * (1 - Omega[portfoliochoice_index])) * K[asset_index] + wagerate * laborincome_states[laborincome_index] * laborincome_shock_states[shockstate_current];

                                                if (cohexo < cohendo[index_ARRLLRRLLP(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)])
                                                {
                                                    save[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = K[0];
                                                    VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = U(cohexo - save[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) + (VFendo[index_ARRLLRRLLP(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - U((cohendo[index_ARRLLRRLLP(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - K[0])));
                                                }

                                                if (cohexo >= cohendo[index_ARRLLRRLLP(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)])
                                                {
                                                    itest = threshold_ii;

                                                    while ((itest < size_asset) && cohexo > cohendo[(index_ARRLLRRLLP(itest, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index))])
                                                    {
                                                        itest++;
                                                    }

                                                    if (itest == size_asset)
                                                    {
                                                        // extrapolation
                                                        vfweight = (cohexo - cohendo[index_ARRLLRRLLP(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (cohendo[index_ARRLLRRLLP(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - cohendo[index_ARRLLRRLLP(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]);
                                                        igridL = size_asset - 2;
                                                        igridH = size_asset - 1;
                                                    }
                                                    else
                                                    {
                                                        // standard interior
                                                        vfweight = (cohexo - cohendo[index_ARRLLRRLLP(itest - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (cohendo[index_ARRLLRRLLP(itest, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - cohendo[index_ARRLLRRLLP(itest - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]);
                                                        igridL = itest - 1;
                                                        igridH = itest - 0;
                                                    }

                                                    VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = inter1d(vfweight, VFendo[index_ARRLLRRLLP(igridL, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VFendo[index_ARRLLRRLLP(igridH, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]);
                                                    save[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = inter1d(vfweight, K[igridL], K[igridH]);

                                                    threshold_ii = min(size_asset - 2, itest);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        temp = 0;
        // std::cout << "rescaling done\n";
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            shockstate_pre = risk_pre_index * size_risk + laborincome_pre_index;

                                            for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
                                            {
                                                tempnext = 0.0;

                                                for (risk_pre_indexnext = 0; risk_pre_indexnext < size_risk; risk_pre_indexnext++)
                                                {
                                                    for (laborincome_pre_indexnext = 0; laborincome_pre_indexnext < size_laborincome; laborincome_pre_indexnext++)
                                                    {
                                                        for (riskshock_pre_indexnext = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                                        {
                                                            for (laborshock_pre_indexnext = 0; laborshock_pre_indexnext < size_shock; laborshock_pre_indexnext++)
                                                            {
                                                                shockstate_current = riskshock_pre_indexnext * size_shock + laborshock_pre_indexnext;
                                                                tempnext += risk_trans[risk_pre_index][risk_pre_indexnext] * laborincome_trans[laborincome_pre_index][laborincome_pre_indexnext] * risk_labor_shock_trans[0][shockstate_current] * VF[index_ARRLLRRLLP(asset_index, risk_pre_indexnext, risk_pre_index, laborincome_pre_indexnext, laborincome_pre_index, riskshock_pre_indexnext, riskshock_pre_index, laborshock_pre_indexnext, laborshock_pre_index, portfoliochoice_index)];
                                                            }
                                                        }
                                                    }
                                                }

                                                // std::cout << tempnext << "\n";
                                                if (portfoliochoice_index == 0)
                                                {
                                                    temp = tempnext;
                                                    itemp = 0;
                                                }

                                                if (tempnext > temp)
                                                {
                                                    temp = tempnext;
                                                    itemp = portfoliochoice_index;
                                                }
                                            }

                                            VF_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, itemp)];
                                            save_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = save[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, itemp)];
                                            Portfolio[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = Omega[itemp];

                                            if (Portfolio[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] > 1)
                                            {
                                                temp += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout << "In while loop// Selection:" << temp << "\n";

        temp = 0;
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {

                                            if (Portfolio[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] > 1)
                                            {
                                                temp += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout << "In while loop// reactivation:" << temp << "\n";

        // std::cout << "port done\n";
        // std::cout << std::setprecision(16) << VF[index_ARRLLRRLLP(5, 5)] << "\n";

        // computing new derivatives and convergence
        critV = 0.0;

        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {

                    for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                    {
                        for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                        {
                            for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                            {
                                for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                {
                                    for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                    {
                                        for (asset_index = 0; asset_index < size_asset; asset_index++)
                                        {

                                            if (asset_index >= 2)
                                            {
                                                dVF_final[index_ARRLLRRLL(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = nderiv(VF_final[index_ARRLLRRLL(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[index_ARRLLRRLL(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                                            }

                                            critV = max(critV, abs(VF_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final_old[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]));

                                            // left corner
                                            dVF_final[index_ARRLLRRLL(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[index_ARRLLRRLL(1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[index_ARRLLRRLL(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[1] - K[0]);
                                            // right corner
                                            dVF_final[index_ARRLLRRLL(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[index_ARRLLRRLL(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[index_ARRLLRRLL(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            for (asset_index = 0; asset_index < size_asset; asset_index++)
                                            {

                                                VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = relaxVF * VF_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] + (1 - relaxVF) * VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            for (asset_index = 0; asset_index < size_asset; asset_index++)
                                            {

                                                if (asset_index >= 2)
                                                {
                                                    dVF[index_ARRLLRRLLP(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = nderiv(VF[index_ARRLLRRLLP(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[index_ARRLLRRLLP(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                                                }

                                                // left corner
                                                dVF[index_ARRLLRRLLP(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[index_ARRLLRRLLP(1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[index_ARRLLRRLLP(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[1] - K[0]);
                                                // right corner
                                                dVF[index_ARRLLRRLLP(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[index_ARRLLRRLLP(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[index_ARRLLRRLLP(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        iter++;
        // std::cout << critV;
        std::cout
            << "iteration=" << iter << ", critV=" << critV << "\n";
    }

    temp = 0;
    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                    {
                        for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                        {
                            for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                            {
                                for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                {
                                    for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                    {

                                        if (Portfolio[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] > 1)
                                        {
                                            temp += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << temp;
}

void SIMULATION(double *save, double *dist, double *capitalout, double K[size_asset])
{
    double *distold, critdist, distverif, weight;

    distold = (double *)calloc((ARRLLRRLL_dim), sizeof(double));
    null(distold, ARRLLRRLL_dim);

    int isave, asset_index, risk_index, risk_indexnext, risk_pre_index, risk_pre_indexnext, laborincome_index, laborincome_indexnext, laborincome_pre_index, laborincome_pre_indexnext, riskshock_index, riskshock_indexnext, riskshock_pre_index, riskshock_pre_indexnext, laborshock_index, laborshock_indexnext, laborshock_pre_index, laborshock_pre_indexnext, iter, shockstate_current_pre, shockstate_next_pre, shockstate_current, shockstate_next;

    critdist = 1.0;
    iter = 0;
    // save[index_ARRLLRRLL(0, 0, 0, 0)] = 0.001;
    while (critdist > epsdist && iter < 100)
    {
        copy(dist, distold, ARRLLRRLL_dim);
        null(dist, ARRLLRRLL_dim);

        // std::cout << critdist << "\n";

        // distribution dynamics
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            if (distold[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] > 0)
                                            {
                                                shockstate_current = risk_index * size_risk + laborincome_index;
                                                shockstate_current_pre = risk_pre_index * size_risk + laborincome_pre_index;

                                                isave = min((int)(floor(getgrid(save[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]))), size_asset - 2);
                                                weight = (save[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - K[isave]) / (K[isave + 1] - K[isave]);
                                                // std::cout << "weight=" << weight << "\n";

                                                for (risk_indexnext = 0; risk_indexnext < size_risk; risk_indexnext++)
                                                {
                                                    for (risk_pre_indexnext = 0; risk_pre_indexnext < size_risk; risk_pre_indexnext++)
                                                    {
                                                        for (laborincome_indexnext = 0; laborincome_indexnext < size_laborincome; laborincome_indexnext++)
                                                        {
                                                            for (laborincome_pre_indexnext = 0; laborincome_pre_indexnext < size_laborincome; laborincome_pre_indexnext++)
                                                            {
                                                                for (riskshock_indexnext = 0; riskshock_indexnext < size_shock; riskshock_indexnext++)
                                                                {
                                                                    for (riskshock_pre_indexnext = 0; riskshock_pre_indexnext < size_shock; riskshock_pre_indexnext++)
                                                                    {
                                                                        for (laborshock_indexnext = 0; laborshock_indexnext < size_shock; laborshock_indexnext++)
                                                                        {
                                                                            for (laborshock_pre_indexnext = 0; laborshock_pre_indexnext < size_shock; laborshock_pre_indexnext++)
                                                                            {
                                                                                shockstate_next_pre = riskshock_pre_indexnext * size_shock + laborshock_pre_indexnext;
                                                                                shockstate_next = riskshock_indexnext * size_shock + laborshock_indexnext;
                                                                                dist[index_ARRLLRRLL(isave, risk_indexnext, risk_pre_indexnext, laborincome_indexnext, laborincome_pre_indexnext, riskshock_indexnext, riskshock_pre_indexnext, laborshock_indexnext, laborshock_pre_indexnext)] += (1.0 - weight) * risk_pre_trans[risk_pre_index][risk_pre_indexnext] * risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * laborincome_pre_trans[laborincome_pre_index][laborincome_pre_indexnext] * risk_labor_shock_pre_trans[shockstate_current_pre][shockstate_next_pre] * risk_labor_shock_trans[shockstate_current][shockstate_next] * distold[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)];
                                                                                dist[index_ARRLLRRLL(min(isave + 1, size_asset - 1), risk_indexnext, risk_pre_indexnext, laborincome_indexnext, laborincome_pre_indexnext, riskshock_indexnext, riskshock_pre_indexnext, laborshock_indexnext, laborshock_pre_indexnext)] += (weight)*risk_pre_trans[risk_pre_index][risk_pre_indexnext] * risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * laborincome_pre_trans[laborincome_pre_index][laborincome_pre_indexnext] * risk_labor_shock_pre_trans[shockstate_current_pre][shockstate_next_pre] * risk_labor_shock_trans[shockstate_current][shockstate_next] * distold[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)];
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
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
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            critdist = (max(critdist, abs(dist[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - distold[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)])));
                                            distverif += dist[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        iter++;
    }
    std::cout << "iteration=" << iter << ", critdist=" << critdist << ", distverify=" << distverif << "\n";

    // *capitalout = 0.0;

    // for (asset_index = 0; asset_index < size_asset; asset_index++)
    // {
    //     for (risk_index = 0; risk_index < size_risk; risk_index++)
    //     {
    //         for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
    //         {
    //             for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    //             {
    //                 *capitalout += dist[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index)] * K[asset_index];
    //             }
    //         }
    //     }
    // }
}

int main()
{
    // MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
    double *VF, *dVF, *save, *cons, *Portfolio;                                // for decision rules
    double capital1, capital0, PIB, critprice, taxL, welfare, rrate, wagerate; // for equilibrium

    std::cout << "Initilization Start\n";
    // Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
    VF = (double *)calloc((ARRLLRRLLP_dim), sizeof(double));   // value function
    dVF = (double *)calloc((ARRLLRRLLP_dim), sizeof(double));  // value function derivative
    save = (double *)calloc((ARRLLRRLLP_dim), sizeof(double)); // value function

    double *VF_final, *dVF_final, *save_final;
    double *distin_final, *distout_final, capitalout; // for simulation

    VF_final = (double *)calloc((ARRLLRRLL_dim), sizeof(double));  // value function
    dVF_final = (double *)calloc((ARRLLRRLL_dim), sizeof(double)); // value function derivative

    save_final = (double *)calloc((ARRLLRRLL_dim), sizeof(double));
    Portfolio = (double *)calloc((ARRLLRRLL_dim), sizeof(double));
    // cons = (double *)calloc((ARRLLRRLLP_dim), sizeof(double));
    distin_final = (double *)calloc((ARRLLRRLL_dim), sizeof(double));
    // distout_final = (double *)calloc((ARRLLRRLL_dim), sizeof(double));

    std::cout << ARRLLRRLLP_dim << "\n";
    null(VF, ARRLLRRLLP_dim);
    null(dVF, ARRLLRRLLP_dim);

    null(VF_final, ARRLLRRLL_dim);
    null(dVF_final, ARRLLRRLL_dim);
    null(save_final, ARRLLRRLL_dim);
    null(Portfolio, ARRLLRRLL_dim);
    null(distin_final, ARRLLRRLL_dim);
    // null(distout_final, ARRLLRRLL_dim);

    int shockstate_current, asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index, tempcount, state_current;

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        K[asset_index] = getlevel(asset_index);
    }

    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    {
        Omega[portfoliochoice_index] = getomega(portfoliochoice_index);
        // std::cout << Omega[portfoliochoice_index] << "\n";
    }

    //  No problem with initilization of Portfolio Choice. All smaller than 1.
    // rrate = 0.040237086402090;
    wagerate = 0.8;
    distin_final[0] = 1.0;
    // taxL=0.3

    // initializing value function and initial derivatives

    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    {
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            shockstate_current = riskshock_index * size_shock + laborshock_index;
                                            VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = U(wagerate * laborincome_states[laborincome_index] * laborincome_shock_states[shockstate_current] + (1 + r_f / 2.0 + (r_f + pi + risk_states[risk_index]) * risk_shock_states[shockstate_current] / 2.0) * K[asset_index]); // REQUIERE TO BE INCREASING IN K (the case here)
                                        }
                                    }
                                }
                            }
                        }
                        // std::cout << VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, portfoliochoice_index)] << "\n";
                    }
                }
            }
        }
    }
    // printf("Policy Computation Start");

    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    {
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                        {
                            for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                            {
                                for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                                {
                                    for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                    {
                                        for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                        {
                                            if (asset_index >= 2)
                                            {
                                                dVF[index_ARRLLRRLLP(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = nderiv(VF[index_ARRLLRRLLP(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[index_ARRLLRRLLP(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                                            }

                                            // left corner
                                            dVF[index_ARRLLRRLLP(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[index_ARRLLRRLLP(1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[index_ARRLLRRLLP(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[1] - K[0]);
                                            // right corner
                                            dVF[index_ARRLLRRLLP(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[index_ARRLLRRLLP(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[index_ARRLLRRLLP(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                                            // std::cout << dVF[index_ARRLLRRLLP(asset_index, risk_index, risk_pre_index, laborincome_index, portfoliochoice_index)] << "\n";
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // printf("Policy Computation Start");

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                    {
                        for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                        {
                            for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                            {
                                for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                {
                                    for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                    {
                                        shockstate_current = risk_index * size_risk + laborincome_index;

                                        VF_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = U(wagerate * laborincome_states[laborincome_index] * laborincome_shock_states[shockstate_current] + (1 + r_f / 2.0 + (r_f + pi + risk_states[risk_index]) * risk_shock_states[shockstate_current] / 2.0) * K[asset_index]); // REQUIERE TO BE INCREASING IN K (the case here)
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // printf("Policy Computation Start");

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                    {
                        for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                        {
                            for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                            {
                                for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                {
                                    for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                    {
                                        if (asset_index >= 2)
                                        {
                                            dVF_final[index_ARRLLRRLL(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = nderiv(VF_final[index_ARRLLRRLL(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[index_ARRLLRRLL(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                                        }

                                        // left corner
                                        dVF_final[index_ARRLLRRLL(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[index_ARRLLRRLL(1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[index_ARRLLRRLL(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[1] - K[0]);
                                        // right corner
                                        dVF_final[index_ARRLLRRLL(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[index_ARRLLRRLL(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[index_ARRLLRRLL(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("Policy Computation Start\n");
    POLICY(VF_final, dVF_final, save_final, VF, dVF, save, Portfolio, K, Omega, wagerate);
    printf("Policy Computation Done\n");
    // SIMULATION(save_final, distin_final, &capital1, K);

    // for (asset_index = 0; asset_index < size_asset; asset_index++)
    // {
    //     std::cout << asset_index << "," << getlevel(asset_index) << ",";

    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
    //     for (risk_index = 0; risk_index < size_risk; risk_index++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLLRRLLP(asset_index, risk_index)]);
    //         std::cout << distout[index_ARRLLRRLLP(asset_index, risk_index)] << ",";
    //     }
    //     std::cout << "\n";
    // }

    // CreateFolder(".\\csv\\");
    // CreateFolder(".\\figure\\");
    std::string common = "18,pe=e-9,std=0.01,premium=" + std::to_string(pi) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_portfoliochoice) + ",rho_c=" + std::to_string(rhopar) + ",Ksize=" + std::to_string(size_asset) + ",Kmax=" + std::to_string(kmax) + ",relaxVF=" + std::to_string(relaxVF) + ",beta=" + std::to_string(betapar) + ",corr=" + std::to_string(corr) + ",Ssize=" + std::to_string(size_risk) + ".csv ";
    std::string filename_dist = "csv\\dist" + common;
    std::string filename_policy = "csv\\policy" + common;
    std::string filename_VF = "csv\\VF" + common;
    std::string filename_Port = "csv\\Portfolio" + common;

    // std::string var = "sometext" + std::to_string(pi);
    // std::cout << var;

    // std::ofstream dfilecsv;
    // dfilecsv.open(filename_dist);
    // dfilecsv << "gridnumber,"
    //          << "capital,";
    // tempcount = 0;

    // for (risk_index = 0; risk_index < size_risk; risk_index++)
    // {
    //     for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
    //     {
    //         for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    //         {
    //             for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
    //             {
    //                 for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
    //                 {
    //                     for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
    //                     {
    //                         for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
    //                         {
    //                             for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
    //                             {
    //                                 dfilecsv << "dist[" << risk_index << " " << risk_pre_index << " " << laborincome_index << laborincome_pre_index << riskshock_index << riskshock_pre_index << laborshock_index << laborshock_pre_index << "]";
    //                                 if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk * size_shock * size_shock * size_shock * size_shock - 1)
    //                                 {
    //                                     dfilecsv << ",";
    //                                 }
    //                                 else
    //                                 {
    //                                     dfilecsv << "\n";
    //                                 }
    //                                 tempcount++;
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // for (asset_index = 0; asset_index < size_asset; asset_index++)
    // {
    //     dfilecsv << asset_index << "," << getlevel(asset_index) << ",";
    //     tempcount = 0;
    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
    //     for (risk_index = 0; risk_index < size_risk; risk_index++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLLRRLLP(asset_index, risk_index)]);
    //         for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
    //         {
    //             for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    //             {
    //                 for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
    //                 {
    //                     for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
    //                     {
    //                         for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
    //                         {
    //                             for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
    //                             {
    //                                 for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
    //                                 {
    //                                     if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk * size_shock * size_shock * size_shock * size_shock - 1)
    //                                     {
    //                                         dfilecsv << distin_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] << ",";
    //                                     }
    //                                     else
    //                                     {
    //                                         dfilecsv << distin_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)];
    //                                     }
    //                                     tempcount++;
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
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

    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                {
                    for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                    {
                        for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                        {
                            for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                            {
                                for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                {
                                    policyfilecsv << "policy[" << risk_index << " " << risk_pre_index << " " << laborincome_index << laborincome_pre_index << riskshock_index << riskshock_pre_index << laborshock_index << laborshock_pre_index << "]";
                                    if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk * size_shock * size_shock * size_shock * size_shock - 1)
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
                        }
                    }
                }
            }
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        policyfilecsv << asset_index << "," << getlevel(asset_index) << ",";
        tempcount = 0;
        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLLRRLLP(asset_index, risk_index)]);
            for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                    {
                        for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                        {
                            for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                            {
                                for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                {
                                    for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                    {
                                        if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk * size_shock * size_shock * size_shock * size_shock - 1)
                                        {
                                            policyfilecsv << save_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] << ",";
                                        }
                                        else
                                        {
                                            policyfilecsv << save_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)];
                                        }
                                        tempcount++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        policyfilecsv << "\n";
    }

    policyfilecsv.close();

    std::ofstream VFfilecsv;
    VFfilecsv.open(filename_VF);
    VFfilecsv << "gridnumber,"
              << "capital,";
    tempcount = 0;

    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                {
                    for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                    {
                        for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                        {
                            for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                            {
                                for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                {
                                    VFfilecsv << "VF[" << risk_index << " " << risk_pre_index << " " << laborincome_index << laborincome_pre_index << riskshock_index << riskshock_pre_index << laborshock_index << laborshock_pre_index << "]";
                                    if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk * size_shock * size_shock * size_shock * size_shock - 1)
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
                        }
                    }
                }
            }
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        VFfilecsv << asset_index << "," << getlevel(asset_index) << ",";
        tempcount = 0;
        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLLRRLLP(asset_index, risk_index)]);
            for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                    {
                        for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                        {
                            for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                            {
                                for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                {
                                    for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                    {
                                        if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk * size_shock * size_shock * size_shock * size_shock - 1)
                                        {
                                            VFfilecsv << VF_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] << ",";
                                        }
                                        else
                                        {
                                            VFfilecsv << VF_final[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)];
                                        }
                                        tempcount++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        VFfilecsv << "\n";
    }

    VFfilecsv.close();

    std::ofstream Portfilecsv;
    Portfilecsv.open(filename_Port);
    Portfilecsv << "gridnumber,"
                << "capital,";
    tempcount = 0;

    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                {
                    for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                    {
                        for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                        {
                            for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                            {
                                for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                {
                                    Portfilecsv << "Port[" << risk_index << " " << risk_pre_index << " " << laborincome_index << laborincome_pre_index << riskshock_index << riskshock_pre_index << laborshock_index << laborshock_pre_index << "]";
                                    if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk * size_shock * size_shock * size_shock * size_shock - 1)
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
                        }
                    }
                }
            }
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        Portfilecsv << asset_index << "," << getlevel(asset_index) << ",";
        tempcount = 0;
        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLLRRLLP(asset_index, risk_index)]);
            for (risk_pre_index = 0; risk_pre_index < size_risk; risk_pre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    for (laborincome_pre_index = 0; laborincome_pre_index < size_laborincome; laborincome_pre_index++)
                    {
                        for (riskshock_index = 0; riskshock_index < size_shock; riskshock_index++)
                        {
                            for (riskshock_pre_index = 0; riskshock_pre_index < size_shock; riskshock_pre_index++)
                            {
                                for (laborshock_index = 0; laborshock_index < size_shock; laborshock_index++)
                                {
                                    for (laborshock_pre_index = 0; laborshock_pre_index < size_shock; laborshock_pre_index++)
                                    {
                                        if (tempcount < size_laborincome * size_laborincome * size_risk * size_risk * size_shock * size_shock * size_shock * size_shock - 1)
                                        {
                                            Portfilecsv << Portfolio[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] << ",";
                                        }
                                        else
                                        {
                                            Portfilecsv << Portfolio[index_ARRLLRRLL(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)];
                                        }
                                        tempcount++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Portfilecsv << "\n";
    }

    Portfilecsv.close();
}
