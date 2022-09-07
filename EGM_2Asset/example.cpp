#include <windows.h>
#include <stdio.h>
#include <iostream>

const int size_j = 100;
#define getomega(x) ((x + 1.0) / (size_j))


int main()
{
    std::cout << getomega(5) << "\n";
}
