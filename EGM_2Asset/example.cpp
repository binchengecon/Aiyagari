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
double pi;
const double rhopar = 3.0;

#define U(x) (pow((x), (1.0 - rhopar)) / (1.0 - rhopar))
// int main()
// {
//     std::string var = "csv\\policy,pi=" + std::to_string(pi) + ".csv";

//     std::ofstream policyfilecsv;
//     policyfilecsv.open(var);
// }

double sstates[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
double strans[7][7] = {
    {0.006262, 0.060934, 0.242398, 0.382924, 0.241063, 0.060262, 0.006157},
    {0.006245, 0.060822, 0.242175, 0.382924, 0.241285, 0.060374, 0.006175},
    {0.006227, 0.060710, 0.241953, 0.382925, 0.241508, 0.060486, 0.006192},
    {0.006210, 0.060597, 0.241730, 0.382925, 0.241730, 0.060597, 0.006210},
    {0.006192, 0.060486, 0.241508, 0.382925, 0.241953, 0.060710, 0.006227},
    {0.006175, 0.060374, 0.241285, 0.382924, 0.242175, 0.060822, 0.006245},
    {0.006157, 0.060262, 0.241063, 0.382924, 0.242398, 0.060934, 0.006262}};

double sstates_nextperiod[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
double strans_nextperiod[7][7];

// int main(int argc, char **argv)
int main()

{

    // int i;
    // pi = 1;
    // std::cout << "You have entered " << argc
    //           << " arguments:"
    //           << "\n";

    // for (i = 0; i < argc; ++i)
    //     std::cout << argv[i] << "\n";

    // std::istringstream iss(argv[1]);
    // iss >> pi;

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

    int i1, i2, i3, i4;
    for (i1 = 0; i1 < 7; i1++)
    {
        for (i2 = 0; i2 < 7; i2++)
        {
            strans_nextperiod[i1][i2] = 0;
        }
    }
}
