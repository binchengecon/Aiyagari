/** -------------------------------------------------
    Alexandre GAILLARD
    endogenous portfolio choice in an aiyagari model
    
    -------
    HEADER ::
    Contain the parameter of interest and key function (utility etc.).
    -------
**/


/** -------------------- OPTION OF THE CODE -------------------- **/

/** ---- MODEL ---- **/
#define OPT_GE              0       // General equilibrium
#define OPT_utility         1       // 1 == CRRA, 2 == EPSTEIN-ZIN
#define OPT_yprocess        1       // 1 == combine pareto + log-normal distribution,  2 == log_normal discretized with Tauchen.
const int share_endo  = 1;

/** PARALLELIZATION **/
#define OPT_OMP             1       // activate parallelization

/** LOG-FILE **/
#define PRINT_SS            2       // 1:: print the results each iteration of policy, 2 :: print the results after the whole loop. 3:: print the result after each iteration of simulation.
#define PRINT_STAT          0       // print the statistics


/** EQUILIBRIUM OPTION **/
double  EQ_QMAX                          = 30.137282;
double  EQ_QMIN                          = 4.135282;
double  EQ_RELAX                         = 0.1;


/** OPTIMIZATION & MATHEMATICAL **/
int MAX_ITER_SIM        = 10000;           // iteration max of the simulation
int MAX_ITER_VAL        = 10000;           // iteration max of the simulation
int INT_NU_compute      = 20;              // compute the share of risky_asset every X iterations (howard improvement)
double EPS_DIST         = 0.0000000001;    // criterion for simulation
double EPS_VALUE        = 0.00000001;      // criterion for value function
double EPS_EQ           = 0.00001;         // criterion for prices.
double EPS_equilibrium  = 0.00001;         // criterion for equilibrium.



#define CRS                 0       // activate simulated method of moments.
#define SOBOLCREATE         0       // (1 = create a sobol, 0 =  create no sobol)
#define nb_para             6
#define nb_moments          5
#define randomnum           1
#define nbsobol             100




/** -------------------- LIBRARIES -------------------- **/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <iostream>
using namespace std;



/** PRICES MODEL -- GUESS **/
double L_tot,GDP,Q_tot,K_tot,R,W,R_guess=0.05,W_guess=2.0,L_guess=0.825338,RQ,Rstar=1.0;



/** -------------------- CALIBRATION AND GRID DEFINITION -------------------- **/

// ---- INVESTMENT ABILITY
#define length_t    1       // type for investment productivity (theta)
double THETA[length_t]              = {0.0};
double P_THETA[length_t][length_t]  = {1.0};
double I_THETA[length_t]            = {1.0};

// ---- INVESTMENT SHOCK ON PRIVATE EQUITY
// Discretization using a Gaussian-Hermite.
#define length_z    3    // define the grid on idiosyncratic productivity shock on investment -- match volatility of PE.
double Z[length_z]   = {0.2254,1.0,1.7746};
double TZ[length_z]  = {1.0/6.0,2.0/3.0,1.0/6.0};
double omega_corr    = 1.0;     // mimic correlation of labor and risky.
double riskfree      = 0.04;
double premium       = 0.10;

/** -------- LABOR INCOME -------- **/
// ---- PERSISTENT (relevant only of we use process 1)
// use a mixture of Pareto + usual component, following Hubmer et al. (2020) [SEE HIS CODE]
#define length_y    7                  // define grid on idiosyncratic productivity shock (productivity + age)
#define length_y2   7                  // number of possible persistent labor income shocks.
const double rho_y  = 0.95;             // in the range of Healthcote et al. (2010) & Hubmer et al. (2020)
const double std_y  = 0.15;             // in the range of Healthcote et al. (2010) & Hubmer et al. (2020)
const double qtop   = 0.9;              // quantile at which the process becomes a pareto.
const double eta    = 3.0;              // eta shape for labor income earnings.
double ebar         = 1.0;
const double Ugrid[length_y] = {0.1,0.25,0.5,0.75,0.9,0.975,0.999};
double Y[length_y], P_Y[length_y][length_y] = {0}, I_Y[length_y];
int P_Y_reduc[length_y][length_y] = {0}, max_P_Y[length_y] = {0};


// ---- TRANSITORY (assume no persistency)
#define length_e    4                   // ability shock on earnings.
const double urate  = 0.075;            // value in Hubmer et al. (2020)
const double u_ben  = 0.0001;           // assumption.
const double std_e  = 0.25;             // in the range of Healthcote et al. (2010) & Hubmer et al. (2020)
double E[length_e], P_E[length_e];




/** -------------------- GRID OF THE STATE SPACE -------------------- **/
#define length_i    1000    // define the grid of cash on hand
#define length_nu   500     // share of risky assets.

// grid space for value
#define length_yt       (length_y*length_t)
#define length_x        (length_i*length_y*length_t)
#define x(i,y,t)        ((length_i*length_y)*(t)+(length_i)*(y)+(i))

// definition of the grid space.
#define echelle  5.0
#define gridspace(i,xmin,xmax,n) (echelle*(exp((log((xmax/echelle)-((xmin/echelle)-1.0))/(n-1))*(i))+((xmin/echelle)-1.0)))
#define invgridspace(x,xmin,xmax,n) (log((x)/echelle - ((xmin/echelle)-1.0))/(log((xmax/echelle)-((xmin/echelle)-1.0))/(n-1.0)))


// GRID DEFINITION
double PSI[length_i]        = {0};
const double PSI_min        =  0.0;
const double PSI_max        =  500;
double A[length_i]          = {0};
const double A_min          =  0.0;
const double A_max          =  PSI_max;
double NUgrid[length_nu]    = {0};
const double NU_min         =  0.0;
const double NU_max         =  1.0;



/** -------------------- PRODUCTION -------------------- **/
const double TFP            = 1.0;                     // TFP level in the final good sector (normalization)
const double deltapar       = 0.06;                    // depreciation rate.
const double alpha          = 0.4;                     // cobb douglas parameter.
double mupar                = 1.0;                     // DRS parameter.
#define production(Qt,Lt)       (pow((Qt),(alpha))*pow((Lt),(1.0-alpha)))           // production function.
#define MPQ(Xt,Lt)              (alpha*mupar*(pow(((Xt)/(Lt)),(alpha-1.0))))        // marginal production function.
#define MPL(Xt,Lt)              ((1.0-alpha)*(pow(((Xt)/(Lt)),(alpha))))            // marginal production function.



/** -------------------- HOUSEHOLD -------------------- **/
double beta                  = 0.90;         // discount factor
const double rhopar          = 1.5;          // CRRA parameter
const double P_die           = 0.04;         // death probability
const double risk_free       = 0.04;         // risk-free interest rate.

// define utility functions:
#define MUc(xx)                     (pow((xx),-rhopar))
#define invMU(u)                    (pow((u),(-(1.0/rhopar))))
#define U(xx)                       ((pow((xx),(1.0-rhopar)))/(1.0-rhopar))
#define invU(uu)                    (pow(((uu)*(1.0-rhopar)),(1.0/(1.0-rhopar))))

// G-CARA specification:
// risk aversion part (function GAMMA)
//double gammapar0    = 1.0;                   // function parameter gamma wealth
//double gammapar1    = 1.5;                   // decreasing parameter for risk aversion.
//double gammapar2    = 0.0;                   // how the elasticity increases with wealth.
//#define alphafun(psi,vartheta)              ( gammapar0*(vartheta)*(1.0/pow(((psi+10)),gammapar1+gammapar2*(psi))) )
//#define GCARA(value,psi,vartheta)           ( (1.0-exp(-alphafun((psi),(vartheta))*(value)))*(1.0/alphafun((psi),(vartheta))) )
//#define invGCARA(gg,psi,vartheta)           ( log(1.0-alphafun((psi),(vartheta))*(gg))/(-alphafun(psi,vartheta)) )

// fixed risk aversion.
#define GCARA(value,psi,vartheta)           ((value))
#define invGCARA(gg,psi,vartheta)           ((gg))



/** -------------------- STATISTICS -------------------- **/


/** -------------------- SMM -------------------- **/

int idum;   int cseed;
char set_moments_name[nb_moments][100], set_params_name[nb_para][10];
char PARA[200]; char SOBOL[200]; char SMM_fileout[200];
double SMM_dist;
double covmat[nb_moments][nb_moments]; // identity matrix.
double set_genmoments[nb_moments];
double set_obsmoments[nb_moments];
double set_params[nb_para];
double set_params_min[nb_para];
double set_params_max[nb_para];
double SOBOL_seq[nbsobol][nb_para];
double PARA_seq[nbsobol][nb_para];




/** -------------------- OUTPUT FILES -------------------- **/

char FILE_log[80];       // to write past try
const char FILE_policy[]      = "OUTPUT/policy.out";
const char FILE_policy2[]     = "OUTPUT/policy2.out";
const char FILE_endopolicy[]  = "OUTPUT/policyendo.out";
const char FILE_distALL[]     = "OUTPUT/distfile.out";
const char FILE_dist[]        = "OUTPUT/distfileALL.out";
const char FILE_savingrate[]  = "OUTPUT/savingrate.out";

// transition files.
const char FILE_transition_checkout[]   = "OUTPUT/TRANSITION/transitioncheckout.out";
const char FILE_transition_welfare[]    = "OUTPUT/TRANSITION/transition_welfare.out";
const char FILE_transition_results[]    = "OUTPUT/TRANSITION/transition_results.out";

