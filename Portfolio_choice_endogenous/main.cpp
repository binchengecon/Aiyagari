/** -------------------------------------------------
    Alexandre GAILLARD - PHILIPP WANGNER  (GIF PAPER, 2020)
    
    -------
    MAIN ::
    This file contains the libraries used and start the code.
        -   main    :: contains the main code
    -------
**/



/** -------------------- FILES AND LIBRAIRIES -------------------- **/
#include "header.hpp"

// MISC FUNCTION //
#include "misc_functions/asa241.hpp"
#include "misc_functions/asa241.cpp"
#include "misc_functions/readinput.cpp"
#include "misc_functions/useful.cpp"
#include "misc_functions/tauchen.cpp"
#include "misc_functions/zbrent_solver.cpp"
#if SOBOLCREATE == 1
#include "calib_functions/sobol.cpp"
#endif

#include "calib_functions/SMMfun.cpp"   // to compute SMM distance
#if CRS == 1
#include "calib_functions/CRS_2.cpp"    // To find best params
#endif
#include "calibration.cpp"
#include "POLICY.cpp"
//#include "SIMULATION.cpp"
#include "steady_state.cpp"
//#include "transition_deterministic.cpp"





int main(int argc, char* argv[]){

/** ---  SET value for observed moments **/
specify_moments();                      // specify a set of moments + name..
specify_params();                       // set a given set of parameters..


SMM_dist = steady_state(set_params);

/** --- ENDOGENEOUS CALIBRATION ACTIVATED **/
#if CRS == 1
    double *Best_PARA;
    Best_PARA = (double *) calloc((nb_para), sizeof(double));

    myCRS(Best_PARA, steady_state);

    printf("CALIBRATION DONE -------- \n \n ");
    for(int i = 0; i < nb_para; i++) {printf("%20.15f\t", Best_PARA[i]);}
#endif   // end of CRS algorithm


printf("\n \n ---- PROGRAM ENDED ---- \n \n"); getchar();

return 0;

}


