#include <windows.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include "tauchen_multi.cpp"
#include "txt/filename.hpp"
double pi;
double rhopar;
const int size_asset = 150; // number of grid points
const int size_portfoliochoice = 100;
const int size_shock = 3;

const int size_laborincome = 3;
const int size_risk = 1;

#define ARRLLP_dim (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * (size_portfoliochoice + 1))
#define ARRLL_dim (size_asset * size_risk * size_risk * size_laborincome * size_laborincome)

#define index_ARRLLP(asset_gridindex, risk_gridindex, risk_pre_gridindex, laborincome_gridindex, laborincome_pre_gridindex, portfoliochoice_gridindex) (((portfoliochoice_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome)) + ((laborincome_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome)) + ((laborincome_gridindex) * (size_asset * size_risk * size_risk)) + ((risk_pre_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))
#define index_ARRLL(asset_gridindex, risk_gridindex, risk_pre_gridindex, laborincome_gridindex, laborincome_pre_gridindex) (((laborincome_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome)) + ((laborincome_gridindex) * (size_asset * size_risk * size_risk)) + ((risk_pre_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))

#define U(x) (pow((x), (1.0 - rhopar)) / (1.0 - rhopar))
// int main()
// {
//     std::string var = "csv\\policy,pi=" + std::to_string(pi) + ".csv";

//     std::ofstream policyfilecsv;
//     policyfilecsv.open(var);
// }

#define ARRLLRRLLP_dim (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock * size_shock * (size_portfoliochoice + 1))
#define ARRLLRRLL_dim (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock * size_shock)

#define index_ARRLLRRLLP(asset_gridindex, risk_gridindex, risk_pre_gridindex, laborincome_gridindex, laborincome_pre_gridindex, riskshock_gridindex, riskshock_pre_gridindex, laborshock_gridindex, laborshock_pre_gridindex, portfoliochoice_gridindex) (((portfoliochoice_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock * size_shock)) + ((laborshock_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock)) + ((laborshock_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock)) + ((riskshock_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock)) + ((riskshock_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome)) + ((laborincome_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome)) + ((laborincome_gridindex) * (size_asset * size_risk * size_risk)) + ((risk_pre_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))
#define index_ARRLLRRLL(asset_gridindex, risk_gridindex, risk_pre_gridindex, laborincome_gridindex, laborincome_pre_gridindex, riskshock_gridindex, riskshock_pre_gridindex, laborshock_gridindex, laborshock_pre_gridindex) (((laborshock_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock * size_shock)) + ((laborshock_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock * size_shock)) + ((riskshock_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome * size_shock)) + ((riskshock_gridindex) * (size_asset * size_risk * size_risk * size_laborincome * size_laborincome)) + ((laborincome_pre_gridindex) * (size_asset * size_risk * size_risk * size_laborincome)) + ((laborincome_gridindex) * (size_asset * size_risk * size_risk)) + ((risk_pre_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))

void null(double *VectorIN, int dim)
{
    int asset_index;
    for (asset_index = 0; asset_index < dim; asset_index++)
    {
        VectorIN[asset_index] = 0;
    }
}

int main(int argc, char **argv)
// int main()

{
    // std::cout << risk_states[0];
    // int i;
    // // pi = 1;
    // std::cout << "You have entered " << argc
    //           << " arguments:"
    //           << "\n";

    // for (i = 0; i < argc; ++i)
    //     std::cout << argv[i] << "\n";

    // std::istringstream iss(argv[1]);
    // iss >> rhopar;

    // std::cout << rhopar <<"\n";
    // std::cout << U(1);

    // int val = std::stoi(argv[1]);
    // double val2 = std::atof(argv[2]);

    // std::cout << pi << "," << val << "," << val2 << "\n";

    // std::string file = std::to_string(val2);
    // std::cout << file;

    // double num = 12345.0000000000001;

    // std::string trimmedString = std::to_string(num).substr(0, std::to_string(num).find(".") + 2 + 1);
    // std::cout << trimmedString << "\n";

    // double youge = 1e16;
    // std::stringstream ss;
    // ss << youge << ".csv";
    // std::string filename = ss.str();                                                    // filename = 1e+16.csv
    // filename.erase(std::remove(filename.begin(), filename.end(), '+'), filename.end()); // removing the '+' sign
    // std::cout << filename;                                                              // 1e16.csv

    // double Omega[2];
    // for (int portfoliochoice_index = 0; portfoliochoice_index < 2; portfoliochoice_index++)
    // {
    //     Omega[portfoliochoice_index] = (sqrt(5.0) + 1.0) / 2.0;
    //     // Omega[portfoliochoice_index] = 0;
    //     // std::cout << Omega[portfoliochoice_index] << "\n";
    //     printf("%.3f\n", Omega[portfoliochoice_index]);
    // }
    // std::cout << U(0.0);

    // int i1, i2, i3, i4;
    // for (i1 = 0; i1 < 7; i1++)
    // {
    //     for (i2 = 0; i2 < 7; i2++)
    //     {
    //         strans_nextperiod[i1][i2] = 0;
    //     }
    // }
    int shockstate_pre, shockstate_current, shockstate_next, asset_index, ii, risk_index, risk_indexnext, risk_pre_index, risk_pre_indexnext, laborincome_index, laborincome_indexnext, laborincome_pre_index, laborincome_pre_indexnext, riskshock_index, riskshock_indexnext, riskshock_pre_index, riskshock_pre_indexnext, laborshock_index, laborshock_indexnext, laborshock_pre_index, laborshock_pre_indexnext, portfoliochoice_index, iter, threshold_ii, Icase, itest, igridL, igridH, itemp;
    double *VF;
    VF = (double *)calloc((ARRLLRRLLP_dim * 100), sizeof(double)); // value function
    // VF = (double *)calloc((ARRLLP_dim), sizeof(double)); // value function
    null(VF, ARRLLP_dim);
    std::cout << ARRLLRRLLP_dim* 100 << std::endl;

    std::cout << index_ARRLLP(0, 0, 0, 0, 0, 0) << std::endl;
    std::cout << VF[index_ARRLLP(0, 0, 0, 0, 0, 11)] << std::endl;
    std::cout << VF[0] << std::endl;
    // VF[index_ARRLLRRLLP(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)];
}
