/*************************************************/
/**          Bin CHENG - 2022                   **/
/**         AIYAGARI MODEL - EGM METHOD         **/
/**                    2022      
/**     Original ContributorL Alexandre GAILLARD **/               
/*************************************************/


// #include "tauchen.cpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <iterator>

#define max(a,b) (((a)>(b))?(a):(b))



// INDEX //
#define maxigrid 500 // define the grid of saving (next period wealth)
#define maxygrid 7

#define ifulldim (maxigrid*maxygrid)
#define inx(igridindex,jclassindex) (((jclassindex)*(maxigrid))+(igridindex))
#define linspace(x0,xmax,n,i) ((i)*(((xmax)-(x0))/(n))+(x0))


// Income process //
const double p_e = 0.6;
const double std_e = 0.4;
const double m_e = 3;
double prod[maxygrid], ytrans[maxygrid][maxygrid], yinv[maxygrid], Labor;


//Convergence criterion
const double epsilon=    0.000001; // Convergence criterion on VFI
//Relaxation parameters
#define relax 0.9
#define relaxsK 0.5
#define relaxsT 0.5



// GRID for assets //
const double Gridmin=0.0;
const double Gridmax=500.0;

const double Echelle1=1.6;
const double grmin=(Gridmin/Echelle1)-1.0;
const double Exponen=log((Gridmax/Echelle1)-grmin)/(maxigrid-1);
#define phi(x) ( Echelle1*(exp(Exponen*(x))+grmin) )
#define phiinv(x) (log((x)/Echelle1-grmin)/Exponen)

double K[maxigrid];



/***************************/
//  CALIBRATION DEFINITION //
/***************************/

// prices //
double wstar;
double rstar;

// Government steady-state //
const double fracG = 0.2; // 20% of gross market G/[f(K,N)+delta*K] = 0.2;
const double taxK = 0.35;
const double taxL_US = 0.36; // implied by the taxK rate to be revenue neutral

// Production parameters //
const double alphapar = 0.36; //production function parameter
const double deltapar = 0.08; //capital depreciation rate

// Preference parameters //
const double rhopar = 3.0; //CRRA parameter
const double betapar = 0.96; //discount factor

// Marginal utilities //
#define MUc(x) (pow((x),-rhopar))
#define inv_MU(u) (pow((u),(-(1/rhopar))))
#define U(x) (pow((x),(1.0-rhopar))/(1.0-rhopar))


double CDFSTDNormal(double x)
{
	//Function CDFSTDNormal: computes the standard normal CDF using Abramowiz and Stegun (1964) approximation


	// constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
	
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
	
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
	
    return 0.5*(1.0 + sign*y);

	
}


//  Note here that:
// The size_t N template parameter is a deduced integral value based upon the array size passed to the template function. 
// Templates parameters can be:
// non-type template parameter;
// type template parameter;
// template template parameter.

template<size_t Nodes>
void tauchenfun(double rho, double m, double amu, double sigma, double grid[Nodes], double Ptransi[Nodes][Nodes])

{

	
	
	int i,j,k;
	
	
	
	//First lets compute unconditiontal variance of yt
	
	double varyt= (sigma*sigma)/(1-rho*rho);
	
	//Compute stddev of yt
	
	double stdyt= sqrt(varyt);
	
	//Define maximum and minimum grid point
	
	double ynodes[Nodes];
	ynodes[Nodes-1]=m*stdyt;
	ynodes[0]=-ynodes[Nodes-1];
	
	//Define interior nodes
	
	double ynodesinterval=(ynodes[Nodes-1]-ynodes[0])/((Nodes-1)*1.0);
	
	for (i=1; i<(Nodes-1); i++) {
		ynodes[i]=ynodes[i-1]+ynodesinterval;
	}
	
	for (i=0; i<Nodes; i++)
	{
		ynodes[i]=ynodes[i]+amu;
	}
	

	for (i=0; i<Nodes; i++) {
		grid[i]=ynodes[i];
	}
	
	
	//Computing transition probability matrix
	
	double transitionMat[Nodes][Nodes];
	
	
	for (j=0; j<Nodes; j++) 
	{
		for (k=1; k<(Nodes-1); k++) 
		{
			
			transitionMat[j][k]=CDFSTDNormal((ynodes[k]-(1-rho)*amu-rho*ynodes[j]+ynodesinterval/2.0)/sigma)-CDFSTDNormal((ynodes[k]-(1-rho)*amu-rho*ynodes[j]-ynodesinterval/2.0)/sigma);
			
		}
		
		transitionMat[j][0]=CDFSTDNormal((ynodes[0]-(1-rho)*amu-rho*ynodes[j]+ynodesinterval/2.0)/sigma);
		transitionMat[j][Nodes-1]=1.0-CDFSTDNormal((ynodes[Nodes-1]-(1-rho)*amu-rho*ynodes[j]-ynodesinterval/2.0)/sigma);
		
	}
	
	
	for (j=0; j<Nodes; j++) 
	{
		for (k=0; k<(Nodes); k++) 
		{
			Ptransi[j][k] = transitionMat[j][k];
		}
	}
    

}

template<size_t dim, size_t dim2>
void inv_distri(double (&invdist)[dim], double (&prob)[dim][dim2])
{
    double tempdist[dim], critdist, sumdist;
    int i, j;
    
    invdist[1] = 1.0;
    
    critdist = 1.0;
    while(critdist > 0.00000001) {
        
        for(i=0;i<dim;i++){
            tempdist[i] = invdist[i];
        }
        
        // compute the invdist //
        for(i = 0; i<dim; i++){
            invdist[i] = 0.0;
        }
        
        for(i = 0; i<dim; i++){
            for(j = 0; j<dim2; j++) {
                invdist[i] += tempdist[j]*prob[j][i];
            }
        }
        
        critdist = 0.0;
        for(i=0;i<dim;i++) {
            critdist = max(abs(invdist[i] - tempdist[i]), critdist);
        }
        
       //printf("%f %f %f %f %f %f %f", invdist[0], invdist[1], invdist[2], tempdist[0], tempdist[1], tempdist[2], critdist); getchar();
        
    }
    
    // renormalize invdist //
    sumdist = 0.0;
    for(i = 0; i<dim; i++) {
        sumdist += invdist[i];
    }
    for(i=0;i<dim;i++) {
        invdist[i] = invdist[i]/sumdist;
    }
    
}

#define GET_VARIABLE_NAME(Variable)  #Variable


template<size_t dim>
void printarr(double (&array)[dim], char* arrayname)
{
    int i;

    for(i=0; i<dim; i++)
    {
        printf("%s[%d]=%f ",arrayname,i, array[i]);
    }
    printf("\n");
}

template<size_t dim, size_t dim2>
void printarr2d(double (&array)[dim][dim2], char* arrayname)
{
    int i, j;

    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim2;j++)
        {
            printf("%s[%d,%d]=%f ",arrayname,i, j, array[i][j]);
        }
        printf("\n");
    }
    printf("\n");

}


void bascule(double *VectorIN, double *VectorOUT, int dim)
{
    int i;
    for(i=0; i<dim; i++){VectorOUT[i]=VectorIN[i];}
}



// MAIN //
int main()
{
    std::cout << "Here we go .\n" ;
    std::cout << "CD " << CDFSTDNormal(0.1) <<"\n";
    int i, y , k  ;


    tauchenfun(p_e, m_e, 0.0, std_e, prod, ytrans);
    inv_distri(yinv, ytrans);

    Labor = 0;
    printf("PRODUCTIVITY SHOCKS\n");

    for(int y=0; y < maxygrid; y++){
        prod[y] = exp(prod[y]);
        Labor = Labor + yinv[y]*prod[y];
    }

    printarr(yinv, GET_VARIABLE_NAME(yinv));
    printarr(prod, GET_VARIABLE_NAME(prod));
    printarr2d(ytrans, GET_VARIABLE_NAME(ytrans));

    // GRID FOR ASSET //
    for(i=0;i<maxigrid;i++)
    {
        K[i]=phi(i);
    }

    
    // MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
    double *VF, *save, *cons;          // for decision rules
    double *distin,*distout;                        // for simulation
    double capital1,capital0,PIB,critprice,taxL,welfare,taxoutL;    // for equilibrium



    // Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
    VF = (double *) calloc((ifulldim), sizeof(double));      // value function
    save = (double *) calloc((ifulldim), sizeof(double));
    cons = (double *) calloc((ifulldim), sizeof(double));
    distin = (double *) calloc((ifulldim), sizeof(double));
    distout = (double *) calloc((ifulldim), sizeof(double));


    // GUESS PRICE & INITIAL DISTRIBUTION //
    rstar = 0.041;          // interest rate
    distin[0] = 1.0;        // initial distribution

    // START BY GUESSING VF //
    // I choose as an initial guess such that current asset level is two times the next asset level //
    for(i = 0; i < maxigrid; i++){
        for(y = 0; y < maxygrid; y++){
            VF[inx(i,y)] = U(prod[y]*wstar + (1+rstar)*K[i])/(1+betapar);                          // REQUIERE TO BE INCREASING IN K (the case here)
            // printf("%f ",VF[inx(i,y)]);
        }
        // printf("\n");
    }
    // printf("\n");


///// START BIG LOOP OVER INTEREST RATE /////
double critere=1.0;          // convergence criterion
int iter = 0;               // iteration 0


// while(critere>epsilon)
{

    // Assume production function of firms is F(K,L) = K^{alpha}*L^{1-alpha}

    wstar = (1.0-alphapar)*pow(alphapar/(rstar+deltapar),alphapar/(1.0-alphapar));

    // SOLVE POLICY FUNCTION //
    // POLICY_EGM(VF,save,cons,taxK,taxL,iter);
    

// INTEGER //
int i,ii,y,ynext,iter,threshold_ii, Icase;


// INITIALIZATION //
double *Kendo, *VFnew, *Kendo_min, *eVF, *deVF, critere, weight, slope1, slope2, tempvf;
VFnew = (double *) calloc((ifulldim), sizeof(double));          // Value function on the next time grid, next iteration
Kendo = (double *) calloc((ifulldim), sizeof(double));       // endogenous grid values
Kendo_min = (double *) calloc((maxygrid), sizeof(double));       // endogenous grid values
eVF = (double *) calloc((ifulldim), sizeof(double));     // expected value function
deVF = (double *) calloc((ifulldim), sizeof(double));    // derivative of the expected value function


bascule(VF,eVF,ifulldim);



// FIRST TIME COMPUTATION OF EXPECTATION + THRESHOLD + ENDO GRID //
for(y = 0; y < maxygrid; y++){

    /** 1. ENDOGENOUS GRID POINT where the borrowing constraint is binding **/
    Kendo_min[y] = (K[0] + inv_MU(deVF[inx(0,y)]) - wstar*prod[y])/(1+rstar); // this is the implied first asset level for which we achieve the borrowing constraint, so all grid point below should achieve the gridpoint (monotonicity of the value function)

 
    for(i = 0; i < maxigrid; i++){

        /** 2. COMPUTE THE IMPLIED CONSUMPTION LEVEL **/
        
        Kendo[inx(i,y)] = (K[i] + inv_MU(deVF[inx(i,y)]) - wstar*prod[y])/(1+rstar);

        if(i==0)
        {
            deVF[inx(i,y)] = (eVF[inx(i+1,y)]-eVF[inx(i,y)])/(K[i+1] - K[i]);

        }

        if(i > 0 && i < (maxigrid-1))
        {
            deVF[inx(i,y)] = deriv(eVF[inx((i-1),y)],eVF[inx(i,y)],eVF[inx((i+1),y)],K[(i-1)],K[i],K[(i+1)]);
        }

        if(i == (maxigrid-1))
        {
            deVF[inx(i,y)] = (eVF[inx((maxigrid-1),y)]-eVF[inx((maxigrid-2),y)])/(K[(maxigrid-1)] - K[(maxigrid-2)]);
        }

    } // end igrid


// Now, we get the endogenous grid! Kendo.
// We just need to finish the interpolation and complete the 


} // end ygrid

while (critere > epsilon && iter < 1000)
{
    /** 3. INTERPOLATE THE VALUE FUNCTION AND COMPUTE CRITERION **/
    for(y = 0; y < maxygrid; y++){
    
        threshold_ii = 0;
        
        for(i = 0; i < maxigrid; i++){
            
            
            if(K[i] < Kendo_min[y]){
                save[inx(i,y)] = K[0];
                tempvf = eVF[inx(0,y)];
            }
            
            
            if(K[i] >= Kendo_min[y]){
            
                ii = max(threshold_ii,0);
                Icase = 0;
                
                while((K[i]>Kendo[inx(ii,y)]) && (ii < maxigrid))
                {
                    if(ii == (maxigrid-1))
                    {
                        Icase = 2;break;
                    }
                    else
                    {
                        ii++;
                    }
                }

                if(Icase == 2){ // case where you extrapolate.
                    slope1 = (eVF[inx((maxigrid-1),y)] - eVF[inx((maxigrid-2),y)])/(Kendo[inx((maxigrid-1),y)] - Kendo[inx((maxigrid-2),y)]);
                    slope2 = (K[(maxigrid-1)] - K[(maxigrid-2)])/(Kendo[inx((maxigrid-1),y)] - Kendo[inx((maxigrid-2),y)]);
 
                    save[inx(i,y)] = (K[i] - Kendo[inx((maxigrid-1),y)])*slope2 + K[(maxigrid-1)];
                    tempvf = (K[i] - Kendo[inx((maxigrid-1),y)])*slope1 + eVF[inx((maxigrid-1),y)];
                }
                
                
                if(Icase == 0){ // normal case
                    weight = (K[i] - Kendo[inx((ii-1),y)])/(Kendo[inx(ii,y)] - Kendo[inx((ii-1),y)]);
                    
                    save[inx(i,y)] = inter1d(weight,K[(ii-1)],K[ii]);
                    tempvf = inter1d(weight,eVF[inx((ii-1),y)],eVF[inx(ii,y)]);
                }
                
                
                // save localisation of K[i] for next grid point //
                threshold_ii = ii; // for next iteration, then set to the previous solution.
            
            }
            
            cons[inx(i,y)] = prod[y]*wstar + K[i]*(1+rstar*) - save[inx(i,y)];
            VFnew[inx(i,y)] = U(cons[inx(i,y)]) + tempvf;
            
            // COMPUTE CRITERION //
            critere = max(critere,fabs(VF[inx(i,y)] - VFnew[inx(i,y)]));
            
            // SAVE THE NEW VFI //
            VF[inx(i,y)] = VFnew[inx(i,y)];




    iter++;

} 

// end equilibrium loop.




    return 0;
}