#include <iostream>
#include <fstream>
#include <string>
#include <direct.h>

int main()
{
    _mkdir("csv");
    std::ofstream p;
    p.open("csv\\output.csv");

    p << "name\n";
    p.close();
}
