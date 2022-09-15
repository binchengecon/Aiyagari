
#include <cmath>
#include <cstdio>
#include "tauchen.cpp"
#include <iostream>

// INDEX //
// #define maxigrid 500 // define the grid of saving (next period wealth)
// #define maxygrid 7

int main()
{
    // TRANSITION + STATE Z OF ENTREPRENEURS //

    const int size_k = 500; // number of grid points

    const double p_e = 0.000001;
    // const double std_e = 0.01;
    const double std_e = 0.015;

    // const double m_e = 3;
    // const int m_e = 6;
    const int m_e = 3;

    const int maxygrid = 2*m_e+1;   // number of productivity classes
    const int size_m = 2*m_e+1;   // number of productivity classes

    double prod[maxygrid], ytrans[maxygrid][maxygrid], yinv[maxygrid], Labor;

    tauchenfun(p_e, m_e, 0.0, std_e, prod, ytrans);
    // inv_distri(yinv, ytrans);

    Labor = 0.0;
    // printf("PRODUCTIVITY SHOCKS\n");
    // for (int y = 0; y < maxygrid; y++)
    // {
    //     printf("%f\t", prod[y]);
    //     prod[y] = exp(prod[y]);
    //     Labor += yinv[y] * prod[y];
    //     // printf("%f\t", prod[y]);
    // }
    // printf("\n");
    // printf("\n");
    // printf("TRANSITION MATRIX FOR PRODUCTIVITY SHOCKS\n");
    // for (int y = 0; y < maxygrid; y++)
    // {
    //     for (int k = 0; k < maxygrid; k++)
    //     {
    //         printf("%f\t", ytrans[y][k]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    printf("For Copy: PRODUCTIVITY SHOCKS\n");
    printf("For Copy: TRANSITION MATRIX FOR PRODUCTIVITY SHOCKS\n");

    // printf("const double nstates[%d] ={", size_m);
    // for (int y = 0; y < maxygrid; y++)
    // {
    //     // printf("%f\t", prod[y]);
    //     // prod[y] = exp(prod[y]);
    //     Labor += yinv[y] * prod[y];

    //     if (y < maxygrid - 1)
    //     {
    //         printf("exp(%f),\t", prod[y]);
    //     }
    //     else
    //     {
    //         printf("exp(%f)", prod[y]);
    //     }
    // }
    // printf("};\n");

    printf("const double risk_states[%d] ={", size_m);
    for (int y = 0; y < maxygrid; y++)
    {
        // printf("%f\t", prod[y]);
        // prod[y] = exp(prod[y]);
        Labor += yinv[y] * prod[y];
        if (y < maxygrid - 1)
        {
            printf("%.8f, ", prod[y]);
        }
        else
        {
            printf("%.8f", prod[y]);
        }
    }
    printf("};\n");

    printf("const double risk_trans[%d][%d] = {\n", size_m, size_m);
    for (int y = 0; y < maxygrid; y++)
    {
        printf("{");
        for (int k = 0; k < maxygrid; k++)
        {
            // printf("%f,\t", ytrans[y][k]);
            if (k < maxygrid - 1)
            {
                printf("%.10f, ", ytrans[y][k]);
            }
            else
            {
                printf("%.10f", ytrans[y][k]);
            }
        }

        if (y < maxygrid - 1)
        {
            printf("},\n");
        }
        else
        {
            printf("}");
        }

        // printf("\n");
    }
    printf("};\n");
}