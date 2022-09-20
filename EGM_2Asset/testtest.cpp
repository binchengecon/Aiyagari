#include "tauchen_multi.cpp"

int main()
{

    int length_z = 3;
    // gaussian - quadrature (weights and nodes)
    gauher(Z_tmp, TZ, length_z);
    sum_PY = 0;
    for (j = 0; j < length_z; j++)
    {
        Z[j] = mean_z + (-pow(std_z, 2.0) / 2.0 + std_z * Z_tmp[j] * sqrt(2.0));
    }
}