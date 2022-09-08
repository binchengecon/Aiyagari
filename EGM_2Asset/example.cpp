#include <windows.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

const double pi = 0.005;
int main()
{
    std::string var = "csv\\policy,pi=" + std::to_string(pi) + ".csv";

    std::ofstream policyfilecsv;
    policyfilecsv.open(var);
}
