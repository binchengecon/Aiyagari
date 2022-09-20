#include "tauchen_multi.cpp"

int main()
{

    // Driver Code
    double AR[2][2] = {{0.000000001, 0}, {0, 0.6}},
           Sigma[2][2] = {{0.2*0.2, 1*0.2*0.16}, {1*0.2*0.16, 0.16*0.16}}, I[4][4];
    double ARAR[4][4], IDARAR[4][4], IDARAR_inv[4][4], Sigma_vec[4], Var_Vec[4];
    double IDARAR_Vec[16], IDARAR_Vec_inv[16];

    Kroneckerproduct(AR, AR, ARAR);
    eye(I);
    MatDeductMat(I, ARAR, IDARAR);
    Vec(IDARAR, IDARAR_Vec);
    InverseMat44(IDARAR_Vec, IDARAR_Vec_inv);
    deVec(IDARAR_Vec_inv, IDARAR_inv);
    Vec2(Sigma, Sigma_vec);
    MatProductVec(IDARAR_inv, Sigma_vec, Var_Vec);

    // deVec(DI, E);

    std::cout << "AR: \n";

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            std::cout << AR[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::cout << "ARAR: \n";

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << ARAR[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::cout << "IDARAR: \n";

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << IDARAR[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::cout << "IDARAR_inv: \n";

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << IDARAR_inv[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Sigma_vec: \n";

    for (int i = 0; i < 4; i++)
    {

        std::cout << Sigma_vec[i] << " ";
    }
    std::cout << "\n";

    std::cout << "Var_Vec: \n";

    for (int i = 0; i < 4; i++)
    {

        std::cout << Var_Vec[i] << " ";
    }
    std::cout << "\n";
    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << C[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << "C2 Done \n";

    // for (int i = 0; i < 16; i++)
    // {
    //     std::cout << D[i] << "\n";
    // }
    // std::cout << "D Done \n";

    // for (int i = 0; i < 4; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         std::cout << E[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    // std::cout << "E Done \n";
}