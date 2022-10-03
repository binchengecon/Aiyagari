// Modify the correlation of risky return and labor income

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
#include "txt/Shock_alex.hpp"
#include "txt/Risk_Labor.hpp"

const int size_asset = 200; // number of grid points
const int size_portfoliochoice = 50;

const double Amin = 0.0;
const double Amax = 500.0;

const double COHmin = 0.0;
const double COHmax = 500.0;

const double betapar = 0.75;
const double alphapar = 0.36;
const double deltapar = 0.08;
const double rhopar = 4;
const double rhopar_w = 3;
const double labor = 1.0219882;

const double epsV = 1.0e-8;
const double epsdist = 1.0e-8;
const double epsK = 1.0e-6;
const double relaxsK = 0.005;
const double relaxVF = 0.000;

// grid constants
const double scale1 = 1.6;
const double grmin = (Amin / scale1) - 1.0;
const double exponen = log((Amax / scale1) - grmin) / (size_asset - 1);

const double pi = 0.03;
const double corr = 1;
const double r_f = 0.04;
const double wagerate = 1.8;

// grid space for value
#define RL_dim (size_laborincome * size_risk)
#define ARL_dim (size_asset * size_laborincome * size_risk)
#define index(i, y, t) ((size_asset * size_laborincome) * (t) + (size_asset) * (y) + (i))

// Utility
#define MUc(x) (pow((x), -rhopar))
#define inv_MU(u) (pow((u), (-(1.0 / rhopar))))
#define U(x) (pow((x), (1.0 - rhopar)) / (1.0 - rhopar))
#define Uw(x) (pow((x + 0.01), (1.0 - rhopar_w)) / (1.0 - rhopar_w))

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

double FOC_Port(double omega, int y, int t, int i, double *G_VF, double A[size_asset], double COH[size_asset], double wagerate)
{
    double G_V_next, weight, COHnext, EG_V = 0.0;
    int icoh, riskshock_index, laborshock_index, shockstate_current, z;

    // constrain individuals to invest if they generate negative wealth!
    // computing the next period certainty equivalent operator.

    for (z = 0; z < length_z; z++)
    {

        // COHnext = (1.0 + (r_f + pi + risk_states[t]) * risk_shock_states[shockstate_current] * omega + r_f * (1 - omega)) * A[i] + wagerate * laborincome_states[y] * laborincome_shock_states[shockstate_current];
        // COHnext = ((1.0 + r_f + pi) * risk_states[t] * risk_shock_states[shockstate_current] * omega + (1 + r_f) * (1 - omega)) * A[i] + wagerate * laborincome_states[y] * laborincome_shock_states[shockstate_current];
        COHnext = ((1.0 + r_f + pi) * risk_states[t] * ZR[z] * omega + (1 + r_f) * (1 - omega)) * A[i] + wagerate * laborincome_states[y] * ZY[z];

        icoh = max(0, min((int)(floor(getgrid(COHnext))), size_asset - 2));
        weight = (COHnext - COH[icoh]) / (COH[icoh + 1] - COH[icoh]);
        if (weight > 1)
        {
            weight = 1.0 + (weight - 1.0) / 1.2;
        } // this is a small correction for potential approx in the extrapolation at the top.
        if (icoh < 0 | weight < 0.0)
        {
            printf("mistake: %d %f psinext = %f", icoh, weight, COHnext);
            getchar();
        }

        // correct for bounds.
        G_V_next = inter1d(weight, G_VF[index(icoh, y, t)], G_VF[index(icoh + 1, y, t)]); // next period VF.
        // EG_V += risk_labor_shock_trans[0][shocstate_current] * G_V_next;
        EG_V += TZ[z] * G_V_next;
    }

    EG_V += Uw(A[i]);
    return (EG_V);
}
void POLICY(double *VF, double *dVF, double *save, double *cons, double *Portfolio, double A[size_asset], double COH[size_asset], double Omega[size_portfoliochoice + 1], double wagerate)
{
    // initialize index
    int iter, yt, y, t, i, tnext, ynext, threshold_ii, ii, ipsi, nu, w;
    double CRIT_VF;
    double EVF_next, dEVF_next, psi_next, weight, V_max, tempvf, wS, wV, max_GVF, Porttemp, val_try, currentcoh;

    // initialize pointer
    double *COH_endo, *EVF_endo, *VF_tilde, *VF_new, *G_VF;
    COH_endo = (double *)calloc((ARL_dim), sizeof(double)); // endogenous grid coh.
    EVF_endo = (double *)calloc((ARL_dim), sizeof(double)); // endogenous grid values function next.
    VF_tilde = (double *)calloc((ARL_dim), sizeof(double)); // value second sub-period.
    VF_new = (double *)calloc((ARL_dim), sizeof(double));   // value new current period.
    G_VF = (double *)calloc((ARL_dim), sizeof(double));     // the certainty equivalent operator.

    // loop for fixed point (contraction mapping theorem)
    CRIT_VF = 100.0;
    iter = 0;

    int stop_transi = 0;

    while ((CRIT_VF > epsV) && (iter < 300) && (stop_transi == 0))
    {

        /** ------------------ STEP 1 ------------------
            -   compute the optimal consumption-saving problem given "asset" and "coh'"
            -------------------------------------------- **/
        // printf("step 1\n");
        for (y = 0; y < size_laborincome; y++)
        {
            for (t = 0; t < size_risk; t++)
            {
                for (i = 0; i < size_asset; i++)
                {
                    // compute expected value next period (to change with bequest later.)
                    EVF_next = 0;
                    dEVF_next = 0;
                    for (ynext = 0; ynext < size_laborincome; ynext++)
                    {
                        for (tnext = 0; tnext < size_risk; tnext++)
                        {
                            EVF_next += betapar * (laborincome_trans[y][ynext] * risk_trans[t][tnext] * (VF[index(i, ynext, tnext)]));
                            dEVF_next += betapar * (laborincome_trans[y][ynext] * risk_trans[t][tnext] * (dVF[index(i, ynext, tnext)]));
                        } // end tnext
                    }     // end ynext

                    // compute the current cash-on-hand.
                    COH_endo[index(i, y, t)] = A[i] + inv_MU(dEVF_next); // unconstrained problem :: psi = a' + c
                    EVF_endo[index(i, y, t)] = U(inv_MU(dEVF_next)) + EVF_next;
                }
            } // end t
        }     // end y

        /** ------------------ STEP 2 ------------------
            -   Rescale the solution.
            -------------------------------------------- **/
        // printf("step 2\n");

        for (y = 0; y < size_laborincome; y++)
        {
            for (t = 0; t < size_risk; t++)
            {
                threshold_ii = 0;

                for (i = 0; i < size_asset; i++)
                {
                    currentcoh = COH[i];

                    // case 1: eat the borrowing limit.
                    if (currentcoh < COH_endo[index(0, y, t)])
                    {
                        save[index(i, y, t)] = A[0];
                        tempvf = EVF_endo[index(0, y, t)];
                    }

                    // case 2: extrapolation.
                    if (currentcoh >= COH_endo[index(size_asset - 1, y, t)])
                    {
                        wV = (EVF_endo[index((size_asset - 1), y, t)] - EVF_endo[index((size_asset - 2), y, t)]) / (COH_endo[index(size_asset - 1, y, t)] - COH_endo[index(size_asset - 2, y, t)]);
                        wS = (A[(size_asset - 1)] - A[(size_asset - 2)]) / (COH_endo[index(size_asset - 1, y, t)] - COH_endo[index(size_asset - 2, y, t)]);

                        save[index(i, y, t)] = (currentcoh - COH_endo[index(size_asset - 1, y, t)]) * wS + A[(size_asset - 1)];
                        tempvf = (currentcoh - COH_endo[index(size_asset - 1, y, t)]) * wV + EVF_endo[index((size_asset - 1), y, t)];
                    }

                    // case 3: interior solution.
                    if (currentcoh < COH_endo[index(size_asset - 1, y, t)] && currentcoh >= COH_endo[index(0, y, t)])
                    {

                        ii = max(threshold_ii, 0);

                        while ((currentcoh > COH_endo[index(ii, y, t)]) && (ii < size_asset - 1))
                        {
                            ii++;
                        }
                        threshold_ii = max(ii - 2, 0);

                        // if (ii == size_asset - 1 && currentcoh > COH_endo[index(ii, y, t)])
                        // {
                        //     printf("mistake 4: %f %f %f %d %d \n", currentcoh, COH_endo[index(size_asset - 2, y, t)], COH_endo[index(size_asset - 1, y, t)], ii, i);
                        //     getchar();
                        // }

                        weight = (currentcoh - COH_endo[index((ii - 1), y, t)]) / (COH_endo[index(ii, y, t)] - COH_endo[index((ii - 1), y, t)]);
                        save[index(i, y, t)] = inter1d(weight, A[(ii - 1)], A[ii]);
                        tempvf = inter1d(weight, EVF_endo[index((ii - 1), y, t)], EVF_endo[index(ii, y, t)]);
                    }

                    cons[index(i, y, t)] = currentcoh - save[index(i, y, t)] - Amin;
                    VF_tilde[index(i, y, t)] = U(cons[index(i, y, t)]) + tempvf;
                    G_VF[index(i, y, t)] = VF_tilde[index(i, y, t)]; // extra value in case we want to make risk aversion type/wealth dep
                }
            } // end t
        }     // end y

        // printf("step 3\n");

        /** ------------------ STEP 3 ------------------
            -   compute the optimal A and share of risky assets.
            -   given the second sub-period solution VF_tilde, cons, save, solve for the optimal k, risky share.
            ---------------------------------------------------------------------------------------------------------- **/
        CRIT_VF = 0;
        for (y = 0; y < size_laborincome; y++)
        {
            for (t = 0; t < size_risk; t++)
            {
                for (i = 0; i < size_asset; i++)
                {
                    // compute the optimal share invested.
                    max_GVF = -100000;
                    for (w = 0; w < size_portfoliochoice + 1; w++)
                    {
                        val_try = FOC_Port(Omega[w], y, t, i, G_VF, A, COH, wagerate);

                        if (val_try > max_GVF)
                        {
                            max_GVF = val_try;
                            Porttemp = Omega[w];
                        }
                    }

                    VF_new[index(i, y, t)] = max_GVF;
                    Portfolio[index(i, y, t)] = Porttemp;

                    // compute the criterion
                    if (i < size_asset)
                    {
                        CRIT_VF = max(CRIT_VF, abs(VF[index(i, y, t)] - VF_new[index(i, y, t)]));
                    }

                    // replace the value (contraction mappin apply)
                    VF[index(i, y, t)] = VF_new[index(i, y, t)];

                } // end i

            } // end t
        }     // end y

        // printf("step 4\n");

        for (y = 0; y < size_laborincome; y++)
        {
            for (t = 0; t < size_risk; t++)
            {
                for (i = 0; i < size_asset; i++)
                {
                    if (i >= 2)
                    {
                        dVF[index(i - 1, y, t)] = nderiv(VF[index(i - 2, y, t)], VF[index(i - 1, y, t)], VF[index(i, y, t)], A[i - 2], A[i - 1], A[i]);
                    }

                    // left corner
                    dVF[index(0, y, t)] = (VF[index(1, y, t)] - VF[index(0, y, t)]) / (A[1] - A[0]);
                    // right corner
                    dVF[index(size_asset - 1, y, t)] = (VF[index(size_asset - 1, y, t)] - VF[index(size_asset - 2, y, t)]) / (A[size_asset - 1] - A[size_asset - 2]);
                }
            }
        }
        iter++;

        printf("ITERATION = %d,  CRIT_VF = %20.15f\n", iter, CRIT_VF); // getchar();

    } // end while

    // free up the pointers.
    free(COH_endo);
    free(VF_new);
    free(VF_tilde);
    free(EVF_endo);
    free(G_VF);
}

void SIMULATION(double *save, double *dist, double *capitalout, double K[size_asset])
{
    double *distold, critdist, distverif, weight;

    distold = (double *)calloc((ARL_dim), sizeof(double));
    null(distold, ARL_dim);

    int isave, asset_index, risk_index, risk_indexnext, laborincome_index, laborincome_indexnext, iter;

    critdist = 1.0;
    iter = 0;

    while (critdist > epsdist && iter < 150)
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
                    if (distold[index(asset_index, risk_index, laborincome_index)] > 0)
                    {

                        isave = min((int)(floor(getgrid(save[index(asset_index, risk_index, laborincome_index)]))), size_asset - 1);
                        weight = (save[index(asset_index, risk_index, laborincome_index)] - K[isave]) / (K[isave + 1] - K[isave]);
                        // std::cout << "weight=" << weight << "\n";

                        for (risk_indexnext = 0; risk_indexnext < size_risk; risk_indexnext++)
                        {
                            for (laborincome_indexnext = 0; laborincome_indexnext < size_laborincome; laborincome_indexnext++)
                            {
                                dist[index(isave, risk_indexnext, laborincome_indexnext)] += (1.0 - weight) * risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * distold[index(asset_index, risk_index, laborincome_index)];
                                dist[index(min(isave + 1, size_asset - 1), risk_indexnext, laborincome_indexnext)] += (weight)*risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * (distold[index(asset_index, risk_index, laborincome_index)]);
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
                    critdist = (max(critdist, abs(dist[index(asset_index, risk_index, laborincome_index)] - distold[index(asset_index, risk_index, laborincome_index)])));
                    distverif += dist[index(asset_index, risk_index, laborincome_index)];
                }
            }
        }

        iter++;
    }

    std::cout << "iteration=" << iter << ", critdist=" << critdist << ", distverify=" << distverif << "\n";

    double error = 0.0;

    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
        {
            error += save[index(0, risk_index, laborincome_index)] * K[asset_index];
        }
    }

    if (error == 0)
    {
        printf("no saving when asset is zero\n ");
    }
}

void printcsv(std::string filename, std::string filetype, double *function)
{

    std::ofstream filecsv;
    filecsv.open(filename);

    int tempcount, laborincome_index, risk_index, asset_index;

    filecsv << "gridnumber,"
            << "capital,";

    tempcount = 0;
    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)

    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            filecsv << filetype << "[" << laborincome_index << " " << risk_index << "]";

            if (tempcount < size_laborincome * size_risk - 1)
            {
                filecsv << ",";
            }
            else
            {
                filecsv << "\n";
            }
            tempcount++;
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        filecsv << asset_index << "," << getlevel(asset_index) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARLP(asset_index, risk_index)]);
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {

                if (risk_index + laborincome_index * size_risk < size_risk * size_laborincome - 1)
                {
                    filecsv << function[index(asset_index, laborincome_index, risk_index)] << ",";
                }

                if (risk_index + laborincome_index * size_risk == size_risk * size_laborincome - 1)
                {
                    filecsv << function[index(asset_index, laborincome_index, risk_index)];
                }
            }
        }

        filecsv << "\n";
    }

    filecsv.close();
}

int main()
{

    double *VF, *dVF, *save, *cons, *Portfolio, *COH, *Omega, *A, *distin, *distout; // for decision rules
    double capital1;                                                                 // for equilibrium

    std::cout << "Initilization Start\n";
    // Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
    VF = (double *)calloc((ARL_dim), sizeof(double));   // value function
    dVF = (double *)calloc((ARL_dim), sizeof(double));  // value function derivative
    save = (double *)calloc((ARL_dim), sizeof(double)); // value function
    cons = (double *)calloc((ARL_dim), sizeof(double)); // value function
    Portfolio = (double *)calloc((ARL_dim), sizeof(double));
    COH = (double *)calloc((ARL_dim), sizeof(double)); // value function
    Omega = (double *)calloc((size_portfoliochoice + 1), sizeof(double));
    A = (double *)calloc((size_asset), sizeof(double));
    distin = (double *)calloc((ARL_dim), sizeof(double));
    distout = (double *)calloc((ARL_dim), sizeof(double));

    // std::cout << ARL_dim << "\n";
    null(VF, ARL_dim);
    null(dVF, ARL_dim);
    null(Portfolio, ARL_dim);
    null(save, ARL_dim);
    null(cons, ARL_dim);

    // null(distout_final, ARL_dim);
    distin[0] = 1.0;

    int shockstate_current, y, i, t, asset_index, risk_index, laborincome_index, portfoliochoice_index, tempcount, state_current;

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        A[asset_index] = getlevel(asset_index);
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        COH[asset_index] = getlevel(asset_index);
    }

    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    {
        Omega[portfoliochoice_index] = getomega(portfoliochoice_index);
        // std::cout << Omega[portfoliochoice_index] << "\n";
    }

    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (asset_index = 0; asset_index < size_asset; asset_index++)
            {

                VF[index(asset_index, risk_index, laborincome_index)] = U(wagerate * laborincome_states[laborincome_index] + (1 + r_f + pi + risk_states[risk_index]) * A[asset_index]); // REQUIERE TO BE INCREASING IN K (the case here)
            }
        }
    }

    for (y = 0; y < size_laborincome; y++)
    {
        for (t = 0; t < size_risk; t++)
        {
            for (i = 0; i < size_asset; i++)
            {
                if (i >= 2)
                {
                    dVF[index(i - 1, y, t)] = nderiv(VF[index(i - 2, y, t)], VF[index(i - 1, y, t)], VF[index(i, y, t)], A[i - 2], A[i - 1], A[i]);
                }

                // left corner
                dVF[index(0, y, t)] = (VF[index(1, y, t)] - VF[index(0, y, t)]) / (A[1] - A[0]);
                // right corner
                dVF[index(size_asset - 1, y, t)] = (VF[index(size_asset - 1, y, t)] - VF[index(size_asset - 2, y, t)]) / (A[size_asset - 1] - A[size_asset - 2]);
            }
        }
    }

    printf("Policy Computation Start\n");
    POLICY(VF, dVF, save, cons, Portfolio, A, COH, Omega, wagerate);
    printf("Policy Computation Done\n");
    printf("Simulation Computation Start\n");
    SIMULATION(save, distin, &capital1, A);
    printf("Simulation Computation Done\n");

    std::string common = "20_change,pi=" + std::to_string(pi) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_portfoliochoice) + ",rho_c=" + std::to_string(rhopar) + ",rho_w=" + std::to_string(rhopar_w) + ",Ksize=" + std::to_string(size_asset) + ",Kmax=" + std::to_string(Amax) + ",relaxVF=" + std::to_string(relaxVF) + ",beta=" + std::to_string(betapar) + ",corr=" + std::to_string(corr) + ",Ssize=" + std::to_string(size_risk) + ".csv ";
    // std::string common = "19.csv ";
    std::string filename_dist = ".\\csv\\dist" + common;
    std::string filename_Port = ".\\csv\\Portfolio" + common;
    std::string filename_policy = ".\\csv\\policy" + common;
    std::string filename_VF = ".\\csv\\VF" + common;

    printcsv(filename_policy, "policy", save);
    printcsv(filename_VF, "VF", VF);
    printcsv(filename_Port, "Port", Portfolio);
    printcsv(filename_dist, "dist", distin);
}