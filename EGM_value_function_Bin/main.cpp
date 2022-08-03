/*************************************************/
/**          Bin CHENG - 2022                   **/
/**         AIYAGARI MODEL - EGM METHOD         **/
/**                    2022
/**     Original ContributorL Alexandre GAILLARD **/
/*************************************************/

#include "header.hpp"
#include "useful.cpp"
#include "tauchen.cpp"
#include "POLICY.cpp"
#include "SIMULATION.cpp"

// MAIN //
int main()
{

    // INITIALIZATION //
    int i, y, k, iter;

    timeval t1, t2;
    double elapsedTime;

    // TRANSITION + STATE Z OF ENTREPRENEURS //
    tauchenfun(p_e, m_e, 0.0, std_e, prod, ytrans);
    inv_distri(yinv, ytrans);

    Labor = 0.0;
    printf("PRODUCTIVITY SHOCKS\n");
    for (int y = 0; y < maxygrid; y++)
    {
        prod[y] = exp(prod[y]);
        Labor += yinv[y] * prod[y];
        printf("%f\t", prod[y]);
    }
    printf("\n");
    printf("\n");
    printf("TRANSITION MATRIX FOR PRODUCTIVITY SHOCKS\n");
    for (int y = 0; y < maxygrid; y++)
    {
        for (int k = 0; k < maxygrid; k++)
        {
            printf("%f\t", ytrans[y][k]);
        }
        printf("\n");
    }
    printf("\n");

    // GRID FOR ASSET //
    for (i = 0; i < maxigrid; i++)
    {
        K[i] = phi(i);
    }

    // MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
    double *VF, *save, *cons;                                          // for decision rules
    double *distin, *distout;                                          // for simulation
    double capital1, capital0, PIB, critprice, taxL, welfare, taxoutL; // for equilibrium

    // Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
    VF = (double *)calloc((ifulldim), sizeof(double)); // value function
    save = (double *)calloc((ifulldim), sizeof(double));
    cons = (double *)calloc((ifulldim), sizeof(double));

    distin = (double *)calloc((ifulldim), sizeof(double));
    distout = (double *)calloc((ifulldim), sizeof(double));

    // GUESS PRICE & INITIAL DISTRIBUTION //
    rstar = 0.041;                                                                              // interest rate
    wstar = (1.0 - alphapar) * pow(alphapar / (rstar + deltapar), alphapar / (1.0 - alphapar)); // wage rate
    distin[0] = 1.0;                                                                            // initial distribution
    taxL = 0.3;                                                                                 // initial taxes

    // START BY GUESSING VF //
    // I choose as an initial guess such that current asset level is two times the next asset level //
    for (i = 0; i < maxigrid; i++)
    {
        for (y = 0; y < maxygrid; y++)
        {
            VF[inx(i, y)] = U(prod[y] * wstar + (1 + rstar) * K[i]) / (1 + betapar); // REQUIERE TO BE INCREASING IN K (the case here)
        }
    }

    /** START EQUILIBRIUM FIXED POINT **/
    // start timer
    gettimeofday(&t1, NULL);

    ///// START BIG LOOP OVER INTEREST RATE /////
    critprice = 1.0; // convergence criterion
    iter = 0;        // iteration 0

    // while (critprice > epsprice)
    // {

    wstar = (1.0 - alphapar) * pow(alphapar / (rstar + deltapar), alphapar / (1.0 - alphapar));

    POLICY_EGM(VF, save, cons, iter);
    SIMULATION(save, distin, distout);
    bascule(distout, distin, ifulldim);

    // } // end equilibrium loop.

    // SAVE IN FILES //
    FILE *pfile, *dfile, *vfile;

    // POLICY FUNCTIONS //
    pfile = fopen(policyfile, "w");
    setbuf(pfile, NULL);
    for (y = 0; y < maxygrid; y++)
    {
        for (i = 0; i < maxigrid; i++)
        {
            fprintf(pfile, "%5d\t%20.15f\t%20.15f\n", i, phi(i), save[inx(i, y)]);
        }
        fprintf(pfile, "\n");
    }
    fclose(pfile);

    // DISTRIBUTION //
    dfile = fopen(distfile, "w");
    setbuf(dfile, NULL);
    for (i = 0; i < maxigrid; i++)
    {
        for (y = 0; y < maxygrid; y++)
        {
            fprintf(dfile, "%5d\t%20.15f\t%20.15f\n", i, phi(i), distout[inx(i, y)]);
        }
        fprintf(dfile, "\n");
    }
    fclose(dfile);

    // VALUE FUNCTION //
    vfile = fopen(valuefile, "w");
    setbuf(vfile, NULL);
    for (i = 0; i < maxigrid; i++)
    {
        for (y = 0; y < maxygrid; y++)
        {
            fprintf(vfile, "%5d\t%20.15f\t%20.15f\n", i, phi(i), VF[inx(i, y)]);
        }
        fprintf(vfile, "\n");
    }
    fclose(vfile);

    _mkdir("csv");

    std::ofstream dfilecsv;
    dfilecsv.open("csv\\dist.csv");
    dfilecsv << "gridnumber,"
             << "capital,"
             << "dist[0],"
             << "dist[1],"
             << "dist[2],"
             << "dist[3],"
             << "dist[4],"
             << "dist[5],"
             << "dist[6],\n";
    for (i = 0; i < maxigrid; i++)
    {
        dfilecsv << i << "," << phi(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < maxygrid; y++)
        { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
            dfilecsv << distout[inx(i, y)] << ",";
        }
        dfilecsv << "\n";
    }

    dfilecsv.close();

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
                  << "policy[6],\n";
    for (i = 0; i < maxigrid; i++)
    {
        policyfilecsv << i << "," << phi(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < maxygrid; y++)
        { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
            policyfilecsv << save[inx(i, y)] << ",";
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
              << "VF[6],\n";
    for (i = 0; i < maxigrid; i++)
    {
        VFfilecsv << i << "," << phi(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < maxygrid; y++)
        { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
            VFfilecsv << VF[inx(i, y)] << ",";
        }
        VFfilecsv << "\n";
    }

    VFfilecsv.close();

    std::ofstream consfilecsv;
    consfilecsv.open("csv\\cons.csv");
    consfilecsv << "gridnumber,"
                << "capital,"
                << "cons[0],"
                << "cons[1],"
                << "cons[2],"
                << "cons[3],"
                << "cons[4],"
                << "cons[5],"
                << "cons[6],\n";
    for (i = 0; i < maxigrid; i++)
    {
        consfilecsv << i << "," << phi(i) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", i, phi(i));
        for (y = 0; y < maxygrid; y++)
        { // fprintf(dfilecsv, "%20.15f,", VF[inx(i, y)]);
            consfilecsv << cons[inx(i, y)] << ",";
        }
        consfilecsv << "\n";
    }

    consfilecsv.close();

    return 0;
}