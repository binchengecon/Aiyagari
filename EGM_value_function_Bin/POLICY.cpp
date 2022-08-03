// SOLVE POLICY FUNCTION //
// POLICY_EGM(VF,save,cons,taxK,taxL,iter);

// POLICY ITERATION //
void POLICY_EGM(double *VF, double *save, double *cons, int iter2)
{

    // INTEGER //
    int i, ii, y, ynext, iter, threshold_ii, Icase;

    // INITIALIZATION //
    double *Kendo, *VFnew, *Kendo_min, *eVF, *deVF, critere, weight, slope1, slope2, tempvf;
    VFnew = (double *)calloc((ifulldim), sizeof(double));     // Value function on the next time grid, next iteration
    Kendo = (double *)calloc((ifulldim), sizeof(double));     // endogenous grid values
    Kendo_min = (double *)calloc((maxygrid), sizeof(double)); // endogenous grid values
    eVF = (double *)calloc((ifulldim), sizeof(double));       // expected value function
    deVF = (double *)calloc((ifulldim), sizeof(double));      // derivative of the expected value function

    // START LOOP OVER DECISION RULES //
    critere = 1.0;
    iter = 0;

    bascule(VF, eVF, ifulldim);

    // FIRST TIME COMPUTATION OF EXPECTATION + THRESHOLD + ENDO GRID //
    for (y = 0; y < maxygrid; y++)
    {
        deVF[inx(0, y)] = (eVF[inx(1, y)] - eVF[inx(0, y)]) / (K[1] - K[0]);

        /** 1. ENDOGENOUS GRID POINT where the borrowing constraint is binding **/
        Kendo_min[y] = (K[0] + inv_MU(deVF[inx(0, y)]) - wstar * prod[y]) / (1 + rstar); // this is the implied first asset level for which we achieve the borrowing constraint, so all grid point below should achieve the gridpoint (monotonicity of the value function)
        for (i = 0; i < maxigrid; i++)
        {
            /** 2. COMPUTE THE IMPLIED CONSUMPTION LEVEL **/

            if (i > 0 && i < (maxigrid - 1))
            {
                deVF[inx(i, y)] = deriv(eVF[inx((i - 1), y)], eVF[inx(i, y)], eVF[inx((i + 1), y)], K[(i - 1)], K[i], K[(i + 1)]);
            }
            if (i == (maxigrid - 1))
            {
                deVF[inx(i, y)] = (eVF[inx((maxigrid - 1), y)] - eVF[inx((maxigrid - 2), y)]) / (K[(maxigrid - 1)] - K[(maxigrid - 2)]);
            }

            Kendo[inx(i, y)] = (K[i] + inv_MU(deVF[inx(i, y)]) - wstar * prod[y]) / (1 + rstar);

        } // end igrid
        // Now, we get the endogenous grid! Kendo.
        // We just need to finish the interpolation and complete the
    } // end ygrid

    while (critere > epsilon && iter < 1000)
    {

        // MAIN LOOP //
        critere = 0.0;

        /** 3. INTERPOLATE THE VALUE FUNCTION AND COMPUTE CRITERION **/
        for (y = 0; y < maxygrid; y++)
        {
            threshold_ii = 0;
            for (i = 0; i < maxigrid; i++)
            {
                if (K[i] < Kendo_min[y])
                {
                    save[inx(i, y)] = K[0];
                    tempvf = eVF[inx(0, y)];
                }
                if (K[i] >= Kendo_min[y])
                {
                    ii = max(threshold_ii, 0);
                    Icase = 0;
                    while ((K[i] > Kendo[inx(ii, y)]) && (ii < maxigrid))
                    {
                        if (ii == (maxigrid - 1))
                        {
                            Icase = 2;
                            break;
                        }
                        else
                        {
                            ii++;
                        }
                    }
                    if (Icase == 2)
                    { // case where you extrapolate.
                        slope1 = (eVF[inx((maxigrid - 1), y)] - eVF[inx((maxigrid - 2), y)]) / (Kendo[inx((maxigrid - 1), y)] - Kendo[inx((maxigrid - 2), y)]);
                        slope2 = (K[(maxigrid - 1)] - K[(maxigrid - 2)]) / (Kendo[inx((maxigrid - 1), y)] - Kendo[inx((maxigrid - 2), y)]);

                        save[inx(i, y)] = (K[i] - Kendo[inx((maxigrid - 1), y)]) * slope2 + K[(maxigrid - 1)];
                        tempvf = (K[i] - Kendo[inx((maxigrid - 1), y)]) * slope1 + eVF[inx((maxigrid - 1), y)];
                    }

                    if (Icase == 0)
                    { // normal case
                        weight = (K[i] - Kendo[inx((ii - 1), y)]) / (Kendo[inx(ii, y)] - Kendo[inx((ii - 1), y)]);

                        save[inx(i, y)] = inter1d(weight, K[(ii - 1)], K[ii]);
                        tempvf = inter1d(weight, eVF[inx((ii - 1), y)], eVF[inx(ii, y)]);
                    }

                    // save localisation of K[i] for next grid point //
                    threshold_ii = ii; // for next iteration, then set to the previous solution.
                }

                cons[inx(i, y)] = prod[y] * wstar + K[i] * (1 + rstar) - save[inx(i, y)];
                VFnew[inx(i, y)] = U(cons[inx(i, y)]) + tempvf;

                // COMPUTE CRITERION //
                critere = max(critere, fabs(VF[inx(i, y)] - VFnew[inx(i, y)]));

                // SAVE THE NEW VFI //
                VF[inx(i, y)] = VFnew[inx(i, y)];

            } // end igrid

        } // end ygrid

        for (y = 0; y < maxygrid; y++)
        {
            for (i = 0; i < maxigrid; i++)
            {

                eVF[inx(i, y)] = 0.0;
                for (ynext = 0; ynext < maxygrid; ynext++)
                {
                    eVF[inx(i, y)] += betapar * ytrans[y][ynext] * VF[inx(i, ynext)];
                }

                if (i >= 2)
                {
                    deVF[inx((i - 1), y)] = deriv(eVF[inx((i - 2), y)], eVF[inx((i - 1), y)], eVF[inx((i), y)], K[(i - 2)], K[(i - 1)], K[(i)]);
                    // deVF[inx((i-1),y)] = (eVF[inx((i),y)] - eVF[inx((i-2),y)])/(K[(i)]-K[(i-2)]);
                    Kendo[inx((i - 1), y)] = (K[(i - 1)] + inv_MU(deVF[inx((i - 1), y)]) - wstar * prod[y]) / (1 + rstar);
                }
            }

            deVF[inx(0, y)] = (eVF[inx(1, y)] - eVF[inx(0, y)]) / (K[1] - K[0]);
            Kendo_min[y] = (K[0] + inv_MU(deVF[inx(0, y)]) - wstar * prod[y]) / (1 + rstar);
            Kendo[inx(0, y)] = (K[0] + inv_MU(deVF[inx(0, y)]) - wstar * prod[y]) / (1 + rstar);

            deVF[inx((maxigrid - 1), y)] = (eVF[inx((maxigrid - 1), y)] - eVF[inx((maxigrid - 2), y)]) / (K[(maxigrid - 1)] - K[(maxigrid - 2)]);
            Kendo[inx((maxigrid - 1), y)] = (K[(maxigrid - 1)] + inv_MU(deVF[inx((maxigrid - 1), y)]) - wstar * prod[y]) / (1 + rstar);
        }

        // printf("CONVERGENCE: %d, %20.15f\n", iter, critere);

        iter++;

        // end of while loop
    }

    free(Kendo_min);
    free(eVF);
    free(deVF);
    free(VFnew);
    free(Kendo);
}