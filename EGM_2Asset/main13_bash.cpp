// add utility by wealth

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <sys/time.h>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <direct.h>
#include <iomanip>
#include <cmath>
#include <limits>
#include <windows.h>
#include <string>
#include <sstream>

const int size_asset = 300; // number of grid points
const int size_risk = 7;    // number of productivity classes
const int size_laborincome = 7;
const int size_portfoliochoice = 100;

#define ARRLP_dim (size_asset * size_risk * size_risk * size_laborincome * (size_portfoliochoice + 1))
#define ARRL_dim (size_asset * size_risk * size_risk * size_laborincome)

#define index_ARRLP(asset_gridindex, risk_gridindex, riskpre_gridindex, laborincome_gridindex, portfoliochoice_gridindex) (((portfoliochoice_gridindex) * (size_asset * size_risk * size_risk * size_laborincome)) + ((laborincome_gridindex) * (size_asset * size_risk * size_risk)) + ((riskpre_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))
#define index_ARRL(asset_gridindex, risk_gridindex, riskpre_gridindex, laborincome_gridindex) (((laborincome_gridindex) * (size_asset * size_risk * size_risk)) + ((riskpre_gridindex) * (size_asset * size_risk)) + ((risk_gridindex) * (size_asset)) + (asset_gridindex))

const double kmin = 0.0;
const double kmax = 200.0;

const double betapar = 0.9;
const double alphapar = 0.36;
const double deltapar = 0.08;
const double rhopar = 3.0;
const double rhopar_w = 3.0;
const double labor = 1.0219882;

const double epsV = 1.0e-8;
const double epsdist = 1.0e-8;
const double epsK = 1.0e-6;
const double relaxsK = 0.005;
const double relaxVF = 0.50;

// grid constants
const double scale1 = 1.6;
const double grmin = (kmin / scale1) - 1.0;
const double exponen = log((kmax / scale1) - grmin) / (size_asset - 1);

const double laborincome_states[7] = {exp(-0.600000000000000), exp(-0.400000000000000), exp(-0.200000000000000), exp(0.000000000000000), exp(0.200000000000000), exp(0.400000000000000), exp(0.600000000000000)};

const double laborincome_trans[7][7] = {
    {0.046746218637144, 0.217937777267117, 0.397822606398702, 0.266386738072197, 0.065169922261456, 0.005754191945237, 0.000182545418147},
    {0.023199661746751, 0.149524091076020, 0.369020347246402, 0.333823905199677, 0.110578117872631, 0.013276146769082, 0.000577730089437},
    {0.010548958644399, 0.093657511915497, 0.312761268311836, 0.382193227897354, 0.171253064028981, 0.027919224002876, 0.001666745199056},
    {0.004387354018187, 0.053538402796357, 0.242163972572887, 0.399820541225137, 0.242163972572887, 0.053538402796357, 0.004387354018187},
    {0.001666745199056, 0.027919224002876, 0.171253064028981, 0.382193227897354, 0.312761268311837, 0.093657511915497, 0.010548958644399},
    {0.000577730089436, 0.013276146769082, 0.110578117872631, 0.333823905199677, 0.369020347246403, 0.149524091076020, 0.023199661746751},
    {0.000182545418147, 0.005754191945237, 0.065169922261456, 0.266386738072197, 0.397822606398702, 0.217937777267117, 0.046746218637144}};

// const double risk_states[7] = {-0.030619, -0.020412, -0.010206, 0.000000, 0.010206, 0.020412, 0.030619};
// const double risk_trans[7][7] = {
//     {0.026240, 0.152924, 0.361483, 0.328567, 0.114742, 0.015266, 0.000778},
//     {0.016044, 0.114742, 0.328567, 0.361483, 0.152924, 0.024700, 0.001539},
//     {0.009452, 0.082835, 0.287445, 0.382789, 0.196114, 0.038437, 0.002929},
//     {0.005362, 0.057531, 0.242024, 0.390166, 0.242024, 0.057531, 0.005362},
//     {0.002929, 0.038437, 0.196114, 0.382789, 0.287445, 0.082835, 0.009452},
//     {0.001539, 0.024700, 0.152924, 0.361483, 0.328567, 0.114742, 0.016044},
//     {0.000778, 0.015266, 0.114742, 0.328567, 0.361483, 0.152924, 0.026240}};

// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.006262, 0.060934, 0.242398, 0.382924, 0.241063, 0.060262, 0.006157},
//     {0.006245, 0.060822, 0.242175, 0.382924, 0.241285, 0.060374, 0.006175},
//     {0.006227, 0.060710, 0.241953, 0.382925, 0.241508, 0.060486, 0.006192},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006192, 0.060486, 0.241508, 0.382925, 0.241953, 0.060710, 0.006227},
//     {0.006175, 0.060374, 0.241285, 0.382924, 0.242175, 0.060822, 0.006245},
//     {0.006157, 0.060262, 0.241063, 0.382924, 0.242398, 0.060934, 0.006262}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.006262, 0.060934, 0.242398, 0.382924, 0.241063, 0.060262, 0.006157},
//     {0.006245, 0.060822, 0.242175, 0.382924, 0.241285, 0.060374, 0.006175},
//     {0.006227, 0.060710, 0.241953, 0.382925, 0.241508, 0.060486, 0.006192},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006192, 0.060486, 0.241508, 0.382925, 0.241953, 0.060710, 0.006227},
//     {0.006175, 0.060374, 0.241285, 0.382924, 0.242175, 0.060822, 0.006245},
//     {0.006157, 0.060262, 0.241063, 0.382924, 0.242398, 0.060934, 0.006262}};

// const double p_e = 0.000001;
// const double p_e = 0.000001;
// const double std_e = 0.01;
// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241730, 0.382925, 0.241730, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241731, 0.382925, 0.241730, 0.060597, 0.006210},
//     {0.006210, 0.060598, 0.241730, 0.382925, 0.241730, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210},
//     {0.006210, 0.060597, 0.241730, 0.382925, 0.241731, 0.060598, 0.006210}};

// const double p_e = 0.0000001;
// const double std_e = 0.01;

// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.0062096854, 0.0605975845, 0.2417303731, 0.3829249254, 0.2417302395, 0.0605975173, 0.0062096749},
//     {0.0062096836, 0.0605975733, 0.2417303508, 0.3829249254, 0.2417302618, 0.0605975285, 0.0062096766},
//     {0.0062096819, 0.0605975621, 0.2417303286, 0.3829249254, 0.2417302840, 0.0605975397, 0.0062096784},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096784, 0.0605975397, 0.2417302840, 0.3829249254, 0.2417303286, 0.0605975621, 0.0062096819},
//     {0.0062096766, 0.0605975285, 0.2417302618, 0.3829249254, 0.2417303508, 0.0605975733, 0.0062096836},
//     {0.0062096749, 0.0605975173, 0.2417302395, 0.3829249254, 0.2417303731, 0.0605975845, 0.0062096854}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.0062096854, 0.0605975845, 0.2417303731, 0.3829249254, 0.2417302395, 0.0605975173, 0.0062096749},
//     {0.0062096836, 0.0605975733, 0.2417303508, 0.3829249254, 0.2417302618, 0.0605975285, 0.0062096766},
//     {0.0062096819, 0.0605975621, 0.2417303286, 0.3829249254, 0.2417302840, 0.0605975397, 0.0062096784},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096784, 0.0605975397, 0.2417302840, 0.3829249254, 0.2417303286, 0.0605975621, 0.0062096819},
//     {0.0062096766, 0.0605975285, 0.2417302618, 0.3829249254, 0.2417303508, 0.0605975733, 0.0062096836},
//     {0.0062096749, 0.0605975173, 0.2417302395, 0.3829249254, 0.2417303731, 0.0605975845, 0.0062096854}};

// const double p_e = 0.00000001;
// const double std_e = 0.01;

// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.0062096806, 0.0605975542, 0.2417303130, 0.3829249254, 0.2417302996, 0.0605975475, 0.0062096796},
//     {0.0062096805, 0.0605975531, 0.2417303108, 0.3829249254, 0.2417303018, 0.0605975486, 0.0062096798},
//     {0.0062096803, 0.0605975520, 0.2417303085, 0.3829249254, 0.2417303041, 0.0605975497, 0.0062096799},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096799, 0.0605975497, 0.2417303041, 0.3829249254, 0.2417303085, 0.0605975520, 0.0062096803},
//     {0.0062096798, 0.0605975486, 0.2417303018, 0.3829249254, 0.2417303108, 0.0605975531, 0.0062096805},
//     {0.0062096796, 0.0605975475, 0.2417302996, 0.3829249254, 0.2417303130, 0.0605975542, 0.0062096806}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.0062096806, 0.0605975542, 0.2417303130, 0.3829249254, 0.2417302996, 0.0605975475, 0.0062096796},
//     {0.0062096805, 0.0605975531, 0.2417303108, 0.3829249254, 0.2417303018, 0.0605975486, 0.0062096798},
//     {0.0062096803, 0.0605975520, 0.2417303085, 0.3829249254, 0.2417303041, 0.0605975497, 0.0062096799},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096799, 0.0605975497, 0.2417303041, 0.3829249254, 0.2417303085, 0.0605975520, 0.0062096803},
//     {0.0062096798, 0.0605975486, 0.2417303018, 0.3829249254, 0.2417303108, 0.0605975531, 0.0062096805},
//     {0.0062096796, 0.0605975475, 0.2417302996, 0.3829249254, 0.2417303130, 0.0605975542, 0.0062096806}};

// const double p_e = 0.000000001;
// const double std_e = 0.01;

// const double risk_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double risk_trans[7][7] = {
//     {0.0062096802, 0.0605975512, 0.2417303070, 0.3829249254, 0.2417303056, 0.0605975505, 0.0062096801},
//     {0.0062096801, 0.0605975511, 0.2417303067, 0.3829249254, 0.2417303059, 0.0605975506, 0.0062096801},
//     {0.0062096801, 0.0605975510, 0.2417303065, 0.3829249254, 0.2417303061, 0.0605975508, 0.0062096801},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096801, 0.0605975508, 0.2417303061, 0.3829249254, 0.2417303065, 0.0605975510, 0.0062096801},
//     {0.0062096801, 0.0605975506, 0.2417303059, 0.3829249254, 0.2417303067, 0.0605975511, 0.0062096801},
//     {0.0062096801, 0.0605975505, 0.2417303056, 0.3829249254, 0.2417303070, 0.0605975512, 0.0062096802}};

// const double riskpre_states[7] = {-0.030000, -0.020000, -0.010000, 0.000000, 0.010000, 0.020000, 0.030000};
// const double riskpre_trans[7][7] = {
//     {0.0062096802, 0.0605975512, 0.2417303070, 0.3829249254, 0.2417303056, 0.0605975505, 0.0062096801},
//     {0.0062096801, 0.0605975511, 0.2417303067, 0.3829249254, 0.2417303059, 0.0605975506, 0.0062096801},
//     {0.0062096801, 0.0605975510, 0.2417303065, 0.3829249254, 0.2417303061, 0.0605975508, 0.0062096801},
//     {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
//     {0.0062096801, 0.0605975508, 0.2417303061, 0.3829249254, 0.2417303065, 0.0605975510, 0.0062096801},
//     {0.0062096801, 0.0605975506, 0.2417303059, 0.3829249254, 0.2417303067, 0.0605975511, 0.0062096801},
//     {0.0062096801, 0.0605975505, 0.2417303056, 0.3829249254, 0.2417303070, 0.0605975512, 0.0062096802}};

const double p_e = 0.000000001;
const double std_e = 0.2;

const double risk_states[7] = {-0.60000000, -0.40000000, -0.20000000, 0.00000000, 0.20000000, 0.40000000, 0.60000000};
const double risk_trans[7][7] = {
    {0.0062096802, 0.0605975512, 0.2417303070, 0.3829249254, 0.2417303056, 0.0605975505, 0.0062096801},
    {0.0062096801, 0.0605975511, 0.2417303067, 0.3829249254, 0.2417303059, 0.0605975506, 0.0062096801},
    {0.0062096801, 0.0605975510, 0.2417303065, 0.3829249254, 0.2417303061, 0.0605975508, 0.0062096801},
    {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
    {0.0062096801, 0.0605975508, 0.2417303061, 0.3829249254, 0.2417303065, 0.0605975510, 0.0062096801},
    {0.0062096801, 0.0605975506, 0.2417303059, 0.3829249254, 0.2417303067, 0.0605975511, 0.0062096801},
    {0.0062096801, 0.0605975505, 0.2417303056, 0.3829249254, 0.2417303070, 0.0605975512, 0.0062096802}};

const double riskpre_states[7] = {-0.60000000, -0.40000000, -0.20000000, 0.00000000, 0.20000000, 0.40000000, 0.60000000};
const double riskpre_trans[7][7] = {
    {0.0062096802, 0.0605975512, 0.2417303070, 0.3829249254, 0.2417303056, 0.0605975505, 0.0062096801},
    {0.0062096801, 0.0605975511, 0.2417303067, 0.3829249254, 0.2417303059, 0.0605975506, 0.0062096801},
    {0.0062096801, 0.0605975510, 0.2417303065, 0.3829249254, 0.2417303061, 0.0605975508, 0.0062096801},
    {0.0062096801, 0.0605975509, 0.2417303063, 0.3829249254, 0.2417303063, 0.0605975509, 0.0062096801},
    {0.0062096801, 0.0605975508, 0.2417303061, 0.3829249254, 0.2417303065, 0.0605975510, 0.0062096801},
    {0.0062096801, 0.0605975506, 0.2417303059, 0.3829249254, 0.2417303067, 0.0605975511, 0.0062096801},
    {0.0062096801, 0.0605975505, 0.2417303056, 0.3829249254, 0.2417303070, 0.0605975512, 0.0062096802}};

const double pi = 0.01;

const double r_f = 0.03;

// Function Definitions:

// Utility
#define MUc(x) (pow((x), -rhopar))
#define inv_MU(u) (pow((u), (-(1.0 / rhopar))))
#define U(x) (pow((x), (1.0 - rhopar)) / (1.0 - rhopar))
#define Uw(x) (pow((x + 0.001), (1.0 - rhopar_w)) / (1.0 - rhopar_w))

// Grid
#define inter1d(x1, y1, y2) ((1.0 - (x1)) * (y1) + (x1) * (y2))
#define getwage(rrate) (1.0 - alphapar) * pow((alphapar / (rrate + deltapar)), (alphapar / (1.0 - alphapar)));
#define getlevel(x) (scale1 * (exp(exponen * (x)) + grmin))
#define getomega(x) ((x + 0.0) / size_portfoliochoice)
#define getgrid(x) (log((x) / scale1 - grmin) / exponen)

// EGM derivatives
#define nderiv(val1, val2, val3, x1, x2, x3) ((1.0 - (x3 - x2) / (x3 - x1)) * ((val3 - val2) / (x3 - x2)) + ((x3 - x2) / (x3 - x1)) * ((val2 - val1) / (x2 - x1)))

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

double K[size_asset], Omega[size_portfoliochoice + 1];

void CreateFolder(const char *path)
{
    if (!CreateDirectory(path, NULL))
    {
        return;
    }
}

void copy(double *VectorIN, double *VectorOUT, int dim)
{
    int asset_index;
    for (asset_index = 0; asset_index < dim; asset_index++)
    {
        VectorOUT[asset_index] = VectorIN[asset_index];
    }
}

void null(double *VectorIN, int dim)
{
    int asset_index;
    for (asset_index = 0; asset_index < dim; asset_index++)
    {
        VectorIN[asset_index] = 0;
    }
}

void POLICY(double *VF_final, double *dVF_final, double *save_final, double *VF, double *dVF, double *save, double *Portfolio, double K[size_asset], double Omega[size_portfoliochoice + 1], double wagerate)
{

    // INITIALIZATION //
    double *Kendo, *VFnew, *Kendo_min, temp, tempnext, dtempnext, *eVF, *deVF, critV, vfweight, slope1, slope2, tempvf, *consendo, *VFendo, *cohendo, cohexo, *VF_final_old;
    VFendo = (double *)calloc((ARRLP_dim), sizeof(double)); // Value function on the next time grid, next iteration
    VFnew = (double *)calloc((ARRLP_dim), sizeof(double));  // Value function on the next time grid, next iteration

    VF_final_old = (double *)calloc((ARRL_dim), sizeof(double)); // Value function on the next time grid, next iteration
    // Kendo = (double *)calloc((ARRLP_dim), sizeof(double));  // endogenous grid values
    // Kendo_min = (double *)calloc((maxygrid), sizeof(double)); // endogenous grid values
    // eVF = (double *)calloc((ARRLP_dim), sizeof(double));       // expected value function
    // deVF = (double *)calloc((ARRLP_dim), sizeof(double));      // derivative of the expected value function
    cohendo = (double *)calloc((ARRLP_dim), sizeof(double));  // Value function on the next time grid, next iteration
    consendo = (double *)calloc((ARRLP_dim), sizeof(double)); // Value function on the next time grid, next iteration

    int asset_index, ii, risk_index, risk_indexnext, riskpre_index, riskpre_indexnext, laborincome_index, laborincome_indexnext, portfoliochoice_index, iter, threshold_ii, Icase, itest, igridL, igridH, itemp;

    iter = 0;

    critV = 10000.0;
    // std::cout << "iter\t"
    //           << "critV\n";

    while (critV > epsV && iter < 250)
    {
        // we need copy to make a separate object

        null(cohendo, ARRLP_dim);
        null(VFendo, ARRLP_dim);

        copy(VF_final, VF_final_old, ARRL_dim);

        // std::cout << std::setprecision(16) << VF[index_ARRLP(5, 5)] << "\n";
        // std::cout << std::setprecision(16) << VFnew[index_ARRLP(5, 5)] << "\n";

        // main EGM computation
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    // try omega here: omega index=portfoliochoice_index
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
                        {

                            tempnext = 0;
                            dtempnext = 0;

                            for (risk_indexnext = 0; risk_indexnext < size_risk; risk_indexnext++)
                            {
                                for (laborincome_indexnext = 0; laborincome_indexnext < size_laborincome; laborincome_indexnext++)
                                {
                                    tempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * VF[index_ARRLP(asset_index, risk_indexnext, risk_index, laborincome_indexnext, portfoliochoice_index)];
                                    dtempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * dVF[index_ARRLP(asset_index, risk_indexnext, risk_index, laborincome_indexnext, portfoliochoice_index)];
                                }
                            }

                            cohendo[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = K[asset_index] + inv_MU(betapar * dtempnext);
                            VFendo[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = U(cohendo[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - K[asset_index]) + Uw(K[asset_index]) + betapar * tempnext;
                        }
                    }
                }
            }
        }
        // std::cout << "EGM done\n";

        // rescalings
        for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
        {

            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        threshold_ii = 0;

                        for (asset_index = 0; asset_index < size_asset; asset_index++)
                        {
                            // method 1: cash on hand
                            cohexo = (1.0 + (r_f + pi + risk_states[risk_index]) * Omega[portfoliochoice_index] + r_f * (1 - Omega[portfoliochoice_index])) * K[asset_index] + wagerate * laborincome_states[laborincome_index];

                            if (cohexo < cohendo[index_ARRLP(0, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)])
                            {
                                save[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = K[0];
                                VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = U(cohexo - save[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]) + Uw(K[asset_index]) + (VFendo[index_ARRLP(0, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - U((cohendo[index_ARRLP(0, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - K[0])));
                            }

                            if (cohexo >= cohendo[index_ARRLP(0, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)])
                            {
                                itest = threshold_ii;

                                while ((itest < size_asset) && cohexo > cohendo[(index_ARRLP(itest, risk_index, riskpre_index, laborincome_index, portfoliochoice_index))])
                                {
                                    itest++;
                                }

                                if (itest == size_asset)
                                {
                                    // extrapolation
                                    vfweight = (cohexo - cohendo[index_ARRLP(size_asset - 2, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]) / (cohendo[index_ARRLP(size_asset - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - cohendo[index_ARRLP(size_asset - 2, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]);
                                    igridL = size_asset - 2;
                                    igridH = size_asset - 1;
                                }
                                else
                                {
                                    // standard interior
                                    vfweight = (cohexo - cohendo[index_ARRLP(itest - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]) / (cohendo[index_ARRLP(itest, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - cohendo[index_ARRLP(itest - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]);
                                    igridL = itest - 1;
                                    igridH = itest - 0;
                                }

                                VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = inter1d(vfweight, VFendo[index_ARRLP(igridL, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)], VFendo[index_ARRLP(igridH, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]);
                                save[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = inter1d(vfweight, K[igridL], K[igridH]);

                                threshold_ii = min(size_asset - 2, itest);
                            }
                        }
                    }
                }
            }
        }
        // std::cout << "rescaling done\n";
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
                        {
                            tempnext = 0.0;
                            for (riskpre_indexnext = 0; riskpre_indexnext < size_risk; riskpre_indexnext++)
                            {
                                tempnext += riskpre_trans[riskpre_index][riskpre_indexnext] * VF[index_ARRLP(asset_index, riskpre_indexnext, riskpre_index, laborincome_index, portfoliochoice_index)];
                            }

                            // std::cout << tempnext << "\n";
                            if (portfoliochoice_index == 0)
                            {
                                temp = tempnext;
                                itemp = 0;
                            }

                            if (tempnext > temp)
                            {
                                temp = tempnext;
                                itemp = portfoliochoice_index;
                            }
                        }

                        VF_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] = VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, itemp)];
                        save_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] = save[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, itemp)];
                        Portfolio[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] = Omega[itemp];
                    }
                }
            }
        }

        // std::cout << "port done\n";
        // std::cout << std::setprecision(16) << VF[index_ARRLP(5, 5)] << "\n";

        // computing new derivatives and convergence
        critV = 0.0;

        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    for (asset_index = 0; asset_index < size_asset; asset_index++)
                    {

                        if (asset_index >= 2)
                        {
                            dVF_final[index_ARRL(asset_index - 1, risk_index, riskpre_index, laborincome_index)] = nderiv(VF_final[index_ARRL(asset_index - 2, risk_index, riskpre_index, laborincome_index)], VF_final[index_ARRL(asset_index - 1, risk_index, riskpre_index, laborincome_index)], VF_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                        }

                        critV = max(critV, abs(VF_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] - VF_final_old[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)]));

                        // left corner
                        dVF_final[index_ARRL(0, risk_index, riskpre_index, laborincome_index)] = (VF_final[index_ARRL(1, risk_index, riskpre_index, laborincome_index)] - VF_final[index_ARRL(0, risk_index, riskpre_index, laborincome_index)]) / (K[1] - K[0]);
                        // right corner
                        dVF_final[index_ARRL(size_asset - 1, risk_index, riskpre_index, laborincome_index)] = (VF_final[index_ARRL(size_asset - 1, risk_index, riskpre_index, laborincome_index)] - VF_final[index_ARRL(size_asset - 2, risk_index, riskpre_index, laborincome_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                    }
                }
            }
        }

        for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (asset_index = 0; asset_index < size_asset; asset_index++)
                        {

                            VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = relaxVF * VF_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] + (1 - relaxVF) * VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)];
                        }
                    }
                }
            }
        }

        for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        for (asset_index = 0; asset_index < size_asset; asset_index++)
                        {

                            if (asset_index >= 2)
                            {
                                dVF[index_ARRLP(asset_index - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = nderiv(VF[index_ARRLP(asset_index - 2, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)], VF[index_ARRLP(asset_index - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)], VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                            }

                            // left corner
                            dVF[index_ARRLP(0, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = (VF[index_ARRLP(1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - VF[index_ARRLP(0, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]) / (K[1] - K[0]);
                            // right corner
                            dVF[index_ARRLP(size_asset - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = (VF[index_ARRLP(size_asset - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - VF[index_ARRLP(size_asset - 2, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                        }
                    }
                }
            }
        }

        iter++;

        std::cout << "iteration=" << iter << ", critV=" << critV << "\n";
    }
}

void SIMULATION(double *save, double *dist, double *capitalout, double K[size_asset])
{
    double *distold, critdist, distverif, weight;

    distold = (double *)calloc((ARRL_dim), sizeof(double));
    null(distold, ARRL_dim);

    int isave, asset_index, risk_index, risk_indexnext, riskpre_index, riskpre_indexnext, laborincome_index, laborincome_indexnext, iter;

    critdist = 1.0;
    iter = 0;
    // save[index_ARRL(0, 0, 0, 0)] = 0.001;
    while (critdist > epsdist && iter < 10000)
    {
        copy(dist, distold, ARRL_dim);
        null(dist, ARRL_dim);

        // std::cout << critdist << "\n";

        // distribution dynamics
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        if (distold[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] > 0)
                        {

                            isave = min((int)(floor(getgrid(save[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)]))), size_asset - 2);
                            weight = (save[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] - K[isave]) / (K[isave + 1] - K[isave]);
                            // std::cout << "weight=" << weight << "\n";

                            for (risk_indexnext = 0; risk_indexnext < size_risk; risk_indexnext++)
                            {
                                for (riskpre_indexnext = 0; riskpre_indexnext < size_risk; riskpre_indexnext++)
                                {
                                    for (laborincome_indexnext = 0; laborincome_indexnext < size_laborincome; laborincome_indexnext++)
                                    {
                                        dist[index_ARRL(isave, risk_indexnext, riskpre_indexnext, laborincome_indexnext)] += (1.0 - weight) * risk_trans[risk_index][risk_indexnext] * riskpre_trans[riskpre_index][riskpre_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * distold[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)];
                                        dist[index_ARRL(min(isave + 1, size_asset - 1), risk_indexnext, riskpre_indexnext, laborincome_indexnext)] += (weight)*risk_trans[risk_index][risk_indexnext] * riskpre_trans[riskpre_index][riskpre_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * (distold[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // convergence
        critdist = 0.0;
        distverif = 0.0;

        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        critdist = (max(critdist, abs(dist[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] - distold[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)])));
                        distverif += dist[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)];
                    }
                }
            }
        }

        iter++;
    }
    std::cout << "iteration=" << iter << ", critdist=" << critdist << ", distverify=" << distverif << "\n";

    *capitalout = 0.0;

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    *capitalout += dist[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] * K[asset_index];
                }
            }
        }
    }
}

int main()
{
    // MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
    double *VF, *dVF, *save, *cons, *Portfolio;                                     // for decision rules
    double capital1, capital0, PIB, critprice, taxL, welfare, rrate, wagerate, coh; // for equilibrium

    // Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
    VF = (double *)calloc((ARRLP_dim), sizeof(double));   // value function
    dVF = (double *)calloc((ARRLP_dim), sizeof(double));  // value function derivative
    save = (double *)calloc((ARRLP_dim), sizeof(double)); // value function

    double *VF_final, *dVF_final, *save_final;
    double *distin_final, *distout_final, capitalout; // for simulation

    VF_final = (double *)calloc((ARRL_dim), sizeof(double));  // value function
    dVF_final = (double *)calloc((ARRL_dim), sizeof(double)); // value function derivative

    save_final = (double *)calloc((ARRL_dim), sizeof(double));
    Portfolio = (double *)calloc((ARRL_dim), sizeof(double));
    // cons = (double *)calloc((ARRLP_dim), sizeof(double));
    distin_final = (double *)calloc((ARRL_dim), sizeof(double));
    distout_final = (double *)calloc((ARRL_dim), sizeof(double));

    null(VF, ARRLP_dim);
    null(dVF, ARRLP_dim);

    null(VF_final, ARRL_dim);
    null(dVF_final, ARRL_dim);
    null(save_final, ARRL_dim);
    null(distin_final, ARRL_dim);
    null(distout_final, ARRL_dim);

    int asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index, tempcount;

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        K[asset_index] = getlevel(asset_index);
    }

    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    {
        Omega[portfoliochoice_index] = getomega(portfoliochoice_index);
        // std::cout << Omega[portfoliochoice_index] << "\n";
    }

    // rrate = 0.040237086402090;
    wagerate = 0.8;
    distin_final[0] = 1.0;
    // taxL=0.3

    // initializing value function and initial derivatives

    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    {
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        coh = wagerate * laborincome_states[laborincome_index] + (1 + r_f + pi + risk_states[risk_index]) * K[asset_index];

                        VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = U(coh / 2.0) + Uw(coh / 2.0); // REQUIERE TO BE INCREASING IN K (the case here)
                        // std::cout << VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] << "\n";
                    }
                }
            }
        }
    }
    // printf("Policy Computation Start");

    for (portfoliochoice_index = 0; portfoliochoice_index < size_portfoliochoice + 1; portfoliochoice_index++)
    {
        for (asset_index = 0; asset_index < size_asset; asset_index++)
        {
            for (risk_index = 0; risk_index < size_risk; risk_index++)
            {
                for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
                {
                    for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                    {
                        if (asset_index >= 2)
                        {
                            dVF[index_ARRLP(asset_index - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = nderiv(VF[index_ARRLP(asset_index - 2, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)], VF[index_ARRLP(asset_index - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)], VF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                        }

                        // left corner
                        dVF[index_ARRLP(0, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = (VF[index_ARRLP(1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - VF[index_ARRLP(0, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]) / (K[1] - K[0]);
                        // right corner
                        dVF[index_ARRLP(size_asset - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] = (VF[index_ARRLP(size_asset - 1, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] - VF[index_ARRLP(size_asset - 2, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                        // std::cout << dVF[index_ARRLP(asset_index, risk_index, riskpre_index, laborincome_index, portfoliochoice_index)] << "\n";
                    }
                }
            }
        }
    }
    // printf("Policy Computation Start");

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    coh = wagerate * laborincome_states[laborincome_index] + (1 + r_f + pi + risk_states[risk_index]) * K[asset_index];
                    VF_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] = U(coh / 2.0) + Uw(coh / 2.0); // REQUIERE TO BE INCREASING IN K (the case here)
                    coh = 0;
                }
            }
        }
    }
    // printf("Policy Computation Start");

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        {
            for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {
                    if (asset_index >= 2)
                    {
                        dVF_final[index_ARRL(asset_index - 1, risk_index, riskpre_index, laborincome_index)] = nderiv(VF_final[index_ARRL(asset_index - 2, risk_index, riskpre_index, laborincome_index)], VF_final[index_ARRL(asset_index - 1, risk_index, riskpre_index, laborincome_index)], VF_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index]);
                    }

                    // left corner
                    dVF_final[index_ARRL(0, risk_index, riskpre_index, laborincome_index)] = (VF_final[index_ARRL(1, risk_index, riskpre_index, laborincome_index)] - VF_final[index_ARRL(0, risk_index, riskpre_index, laborincome_index)]) / (K[1] - K[0]);
                    // right corner
                    dVF_final[index_ARRL(size_asset - 1, risk_index, riskpre_index, laborincome_index)] = (VF_final[index_ARRL(size_asset - 1, risk_index, riskpre_index, laborincome_index)] - VF_final[index_ARRL(size_asset - 2, risk_index, riskpre_index, laborincome_index)]) / (K[size_asset - 1] - K[size_asset - 2]);
                }
            }
        }
    }
    printf("Policy Computation Start\n");
    POLICY(VF_final, dVF_final, save_final, VF, dVF, save, Portfolio, K, Omega, wagerate);
    printf("Policy Computation Done\n");
    SIMULATION(save_final, distin_final, &capital1, K);

    // for (asset_index = 0; asset_index < size_asset; asset_index++)
    // {
    //     std::cout << asset_index << "," << getlevel(asset_index) << ",";

    //     //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
    //     for (risk_index = 0; risk_index < size_risk; risk_index++)
    //     { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLP(asset_index, risk_index)]);
    //         std::cout << distout[index_ARRLP(asset_index, risk_index)] << ",";
    //     }
    //     std::cout << "\n";
    // }

    CreateFolder(".\\csv\\");
    CreateFolder(".\\figure\\");

    std::string filename_common = "13,pe=e-9,std=0.2,premium=" + std::to_string(pi) + ",wage=" + std::to_string(wagerate) + ",rf=" + std::to_string(r_f) + ",Psize=" + std::to_string(size_portfoliochoice) + ",rho_c=" + std::to_string(rhopar) + ",rho_w=" + std::to_string(rhopar_w) + ",Ksize=" + std::to_string(size_asset) + ",relaxVF=" + std::to_string(relaxVF) + ",beta=" + std::to_string(betapar) + ".csv ";
    std::string filename_dist = "csv\\dist" + filename_common;
    std::string filename_policy = "csv\\policy" + filename_common;
    std::string filename_VF = "csv\\VF" + filename_common;
    std::string filename_Port = "csv\\Portfolio" + filename_common;

    // std::string var = "sometext" + std::to_string(pi);
    // std::cout << var;

    std::ofstream dfilecsv;
    dfilecsv.open(filename_dist);
    dfilecsv << "gridnumber,"
             << "capital,";
    tempcount = 0;

    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                dfilecsv << "dist[" << risk_index << " " << riskpre_index << " " << laborincome_index << "]";
                if (tempcount < size_laborincome * size_risk * size_risk - 1)
                {
                    dfilecsv << ",";
                }
                else
                {
                    dfilecsv << "\n";
                }
                tempcount++;
            }
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        dfilecsv << asset_index << "," << getlevel(asset_index) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLP(asset_index, risk_index)]);
            for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {

                    if (risk_index + riskpre_index * size_risk + laborincome_index * size_risk * size_laborincome < size_risk * size_risk * size_laborincome - 1)
                    {
                        dfilecsv << distin_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] << ",";
                    }
                    else
                    {
                        dfilecsv << distin_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)];
                    }
                }
            }
        }

        dfilecsv << "\n";
    }

    dfilecsv.close();

    std::ofstream policyfilecsv;
    policyfilecsv.open(filename_policy);
    policyfilecsv << "gridnumber,"
                  << "capital,";
    tempcount = 0;

    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                policyfilecsv << "policy[" << risk_index << " " << riskpre_index << " " << laborincome_index << "]";
                if (tempcount < size_laborincome * size_risk * size_risk - 1)
                {
                    policyfilecsv << ",";
                }
                else
                {
                    policyfilecsv << "\n";
                }
                tempcount++;
            }
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        policyfilecsv << asset_index << "," << getlevel(asset_index) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLP(asset_index, risk_index)]);
            for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {

                    if (risk_index + riskpre_index * size_risk + laborincome_index * size_risk * size_laborincome < size_risk * size_risk * size_laborincome - 1)
                    {
                        policyfilecsv << save_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] << ",";
                    }
                    else
                    {
                        policyfilecsv << save_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)];
                    }
                }
            }
        }

        policyfilecsv << "\n";
    }

    policyfilecsv.close();

    std::ofstream VFfilecsv;
    VFfilecsv.open(filename_VF);
    VFfilecsv << "gridnumber,"
              << "capital,";
    tempcount = 0;

    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                VFfilecsv << "VF[" << risk_index << " " << riskpre_index << " " << laborincome_index << "]";
                if (tempcount < size_laborincome * size_risk * size_risk - 1)
                {
                    VFfilecsv << ",";
                }
                else
                {
                    VFfilecsv << "\n";
                }
                tempcount++;
            }
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        VFfilecsv << asset_index << "," << getlevel(asset_index) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLP(asset_index, risk_index)]);
            for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {

                    if (risk_index + riskpre_index * size_risk + laborincome_index * size_risk * size_laborincome < size_risk * size_risk * size_laborincome - 1)
                    {
                        VFfilecsv << VF_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] << ",";
                    }
                    else
                    {
                        VFfilecsv << VF_final[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)];
                    }
                }
            }
        }

        VFfilecsv << "\n";
    }

    VFfilecsv.close();

    std::ofstream Portfilecsv;
    Portfilecsv.open(filename_Port);
    Portfilecsv << "gridnumber,"
                << "capital,";
    tempcount = 0;

    for (risk_index = 0; risk_index < size_risk; risk_index++)
    {
        for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
        {
            for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
            {
                Portfilecsv << "Portfolio[" << risk_index << " " << riskpre_index << " " << laborincome_index << "]";
                if (tempcount < size_laborincome * size_risk * size_risk - 1)
                {
                    Portfilecsv << ",";
                }
                else
                {
                    Portfilecsv << "\n";
                }
                tempcount++;
            }
        }
    }

    for (asset_index = 0; asset_index < size_asset; asset_index++)
    {
        Portfilecsv << asset_index << "," << getlevel(asset_index) << ",";

        //  fprintf(dfilecsv, "%5d\t,%20.15f\t,", asset_index, phi(asset_index));
        for (risk_index = 0; risk_index < size_risk; risk_index++)
        { // fprintf(dfilecsv, "%20.15f,", VF[index_ARRLP(asset_index, risk_index)]);
            for (riskpre_index = 0; riskpre_index < size_risk; riskpre_index++)
            {
                for (laborincome_index = 0; laborincome_index < size_laborincome; laborincome_index++)
                {

                    if (risk_index + riskpre_index * size_risk + laborincome_index * size_risk * size_laborincome < size_risk * size_risk * size_laborincome - 1)
                    {
                        Portfilecsv << Portfolio[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)] << ",";
                    }
                    else
                    {
                        Portfilecsv << Portfolio[index_ARRL(asset_index, risk_index, riskpre_index, laborincome_index)];
                    }
                }
            }
        }

        Portfilecsv << "\n";
    }

    Portfilecsv.close();
}
