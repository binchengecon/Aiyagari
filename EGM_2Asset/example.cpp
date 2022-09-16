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

// int main()
// {
//     std::string var = "csv\\policy,pi=" + std::to_string(pi) + ".csv";

//     std::ofstream policyfilecsv;
//     policyfilecsv.open(var);
// }

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
    double Omega[2];
    for (int portfoliochoice_index = 0; portfoliochoice_index < 2; portfoliochoice_index++)
    {
        Omega[portfoliochoice_index] = (sqrt(5.0) + 1.0) / 2.0;
        // Omega[portfoliochoice_index] = 0;
        // std::cout << Omega[portfoliochoice_index] << "\n";
        printf("%.3f\n", Omega[portfoliochoice_index]);
    }
}
