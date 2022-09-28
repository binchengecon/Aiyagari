/** -------------------------------------------------
    by Alexandre GAILLARD (2021)  --  GIF PAPER
    
    -------
    Numerical Functions ::
        - set of various function needed to compute the processes of interest.
            - normalinv, normalpdf, normalcdf
            - Qpareto
            - Ztrans (to convert a process into a double process)
            - Gaussian Hermite Quadrature (gauher) translated from Numerical Recipe.
    -------
**/

// useful numbers:
const double PI       = 3.141592653589793238462643383279502884197;
const double PIO2     = 1.57079632679489661923132169163975144209858;
const double TWOPI    = 6.283185307179586476925286766559005768394;
const double SQRT2    = 1.41421356237309504880168872420969807856967;
const double EULER    = 0.5772156649015328606065120900824024310422;


/** MATHEMATICAL FUNCTIONS **/
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define interpol(x,y,z) (y+(x-floor(x))*(z-y))
#define inter1d(x1,y1,y2) ((1.0-(x1))*(y1)+(x1)*(y2))
#define inter2d(x1,x2,y11,y21,y12,y22) ((1.0-(x2))*((1.0-(x1))*(y11)+(x1)*(y21))+(x2)*((1.0-(x1))*(y12)+(x1)*(y22)))
double CosineInterpolate(double mu, double y1,double y2){
   double mu2;

   mu2 = (1-cos(mu*PI))/2;
   return(y1*(1-mu2)+y2*mu2);
}


/** GRID SPACING **/
#define linspace(x0,xmax,n,i) ((i)*(((xmax)-(x0))/(n-1))+(x0)) // transform a grid into a value
#define invlinspace(x0,xmax,n,x) ((((x)-(x0))/((xmax)-(x0)))*(n-1)) // transform a value into a grid
#define expspace(i,xmin,xmax,echelle,n) ( echelle*(exp((log((xmax/echelle)-((xmin/echelle)-1.0))/(n-1))*(i))+((xmin/echelle)-1.0)) ) // transform a grid into a value
#define invexpspace(x,xmin,xmax,echelle,n) ( log((x)/echelle - ((xmin/echelle)-1.0))/(log((xmax/echelle)-((xmin/echelle)-1.0))/(n-1)) ) // transform a value into a grid
// this spacing mix an exponential and then an equispaced grid //
#define fungrid(i,xmin1,xmax1,xmin2,xmax2,echelle,n1,n2) ((i)<(n1))?(echelle*(exp((log((xmax1/echelle)-((xmin1/echelle)-1.0))/(n1-1))*(i))+((xmin1/echelle)-1.0))):((i-n1+1)*(((xmax2)-(xmin2))/(n2))+(xmin2))
#define invfungrid(x,xmin1,xmax1,xmin2,xmax2,echelle,n1,n2) ((x)<=(xmax1))?(log((x)/echelle - ((xmin1/echelle)-1.0))/(log((xmax1/echelle)-((xmin1/echelle)-1.0))/(n1-1))):((((x)-(xmin2))/((xmax2)-(xmin2)))*(n2)+n1-1)


/** COMPUTE DERIVATIVE **/
#define deriv_corner(x1,x2,val1,val2) ((val2 - val1)/(x2-x1))
#define deriv(x1,x2,x3,val1,val2,val3) ((1.0 - (x3 - x2)/(x3 - x1))*((val3 - val2)/(x3-x2)) + ((x3 - x2)/(x3 - x1))*((val2 - val1)/(x2-x1)))



/** SIGN FUNCTION **/
template<class T>
inline const T SIGN(const T &a, const T &b)
    {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
    {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
    {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}


/** TIMER FUNCTION **/
double timer_fun(timeval t1, timeval t2){
    double elapsedTime;
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
    
    return (elapsedTime/1000);
}


double CubicInterpolate(double mu,double y0,double y1,double y2,double y3){
   double a0,a1,a2,a3,mu2;

   mu2 = mu*mu;
   a0 = y3 - y2 - y0 + y1;
   a1 = y0 - y1 - a0;
   a2 = y2 - y0;
   a3 = y1;

   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}



double HermiteInterpolate(double mu,double y0,double y1,double y2,double y3,double tension,double bias){
   double m0,m1,mu2,mu3;
   double a0,a1,a2,a3;

	mu2 = mu * mu;
	mu3 = mu2 * mu;
   m0  = (y1-y0)*(1+bias)*(1-tension)/2;
   m0 += (y2-y1)*(1-bias)*(1-tension)/2;
   m1  = (y2-y1)*(1+bias)*(1-tension)/2;
   m1 += (y3-y2)*(1-bias)*(1-tension)/2;
   a0 =  2*mu3 - 3*mu2 + 1;
   a1 =    mu3 - 2*mu2 + mu;
   a2 =    mu3 -   mu2;
   a3 = -2*mu3 + 3*mu2;

   return(a0*y1+a1*m0+a2*m1+a3*y2);
}


/** QUADRATIC INTERPOLATION **/
double interQuad1d(const double dx, const double h0, const double h1, const double h2) // dx lies in the interval [i-0.5; i+0.5] (distance is equal to 1), we fake f(i-0.5);
{
double h05, h15, a, b, c, fapprox;

h05 = (h1+h0)/2.0;
h15 = (h2+h1)/2.0;

c = h05;
b = 2.0*(h1-h05);
a = h15-2.0*h1+h05;

fapprox = a*dx*dx+b*dx+c;

return fapprox;

}



/** TRILINEAR INTERPOLATION **/
double inter3d(const double dx, const double dy, const double dz, const double c000, const double c001, const double c010, const double c011, const double c100, const double c101, const double c110, const double c111)
{
    double c00, c01, c10, c11, c0, c1, c;
    
    c00 = c000*(1.0-dx)+c100*dx;
    c01 = c001*(1.0-dx)+c101*dx;
    c10 = c010*(1.0-dx)+c110*dx;
    c11 = c011*(1.0-dx)+c111*dx;
    
    c0 = c00*(1.0-dy)+c10*dy;
    c1 = c01*(1.0-dy)+c11*dy;
    
    c = c0*(1.0-dz)+c1*dz;
    
    return c;
}



void makezero(double *VectorIN, int dim)
{
    int i;
    for(i=0; i<dim; i++){VectorIN[i]=0.0;}
}



/**  FIND THE WEIGHT FOR INTERPOLATION **/
double weightinter(double x, double xgrid, int *ixgrid, double vector[], int dimweight)
{
    double dxgrid;
    
    // PUT THESE VALUE ON A GRID //
    *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
    *ixgrid=max(0,*ixgrid);
    
    if (*ixgrid>=(dimweight-1))
    {
        xgrid=(dimweight-1)-0.000000001;
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    if (*ixgrid<=0)
    {
        xgrid=0.000000001;
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    dxgrid=(x-vector[*ixgrid])/(vector[(*ixgrid+1)]-vector[*ixgrid]);
    
    return dxgrid;
    
}




/** COMPARE VECTORS **/
int comparefun2(const void* a, const void* b)
{
    double* da = (double*)a;
    double* db = (double*)b;
    int diff1 = (da[0] > db[0]) - (da[0] < db[0]);
    if (diff1 != 0) return diff1;
    return (da[1] > db[1]) - (da[1] < db[1]);
}




/** SHIFTER **/
void shft2(double &a, double &b, const double c){a=b;b=c;}
void shft3(double &a, double &b, double &c, const double d){a=b;b=c;c=d;}



/** BASCULE FUNCTIONS **/
void bascule(double *VectorIN, double *VectorOUT, int dim)
{
    int i;
    for(i=0; i<dim; i++){VectorOUT[i]=VectorIN[i];}
}

void bascule_zero(double *VectorIN, double *VectorOUT, int dim)
{
    int i;
    for(i=0; i<dim; i++){
        VectorOUT[i]=VectorIN[i];
        VectorIN[i] = 0.0;
    }
}



/** GET RESIDUALS OF THE HISTOGRAM **/
void weighthist(double x, double xgrid, double *residout, int *ixgrid, double vector[], int dimweight)
{
    
    // PUT THESE VALUE ON A GRID //
    *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
    *ixgrid=max(0,*ixgrid);
    
    if (*ixgrid>=(dimweight-1))
    {
        printf("getresid: (ixgrid>(dimweight-1))");getchar();
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    if (*ixgrid<=0)
    {
        printf("getresid: (ixgrid<0)");getchar();
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    *residout=(x-vector[*ixgrid])/(vector[(*ixgrid+1)]-vector[*ixgrid]);
}



void gethist(double x, double *weight, int *index, double vector[], int dim, double valmax){
    double xval;
    xval    = min(x,valmax-0.0000001);
    *index  = (int)(floor(invexpspace(xval,0,valmax,1.0,dim)));
    *weight = (xval - vector[*index])/(vector[(*index+1)] - vector[*index]);
    if(*index<0){printf("getresid: (index<0), %f %f %f and %d", xval, x, *weight, *index);getchar();}
    if(*index>=(dim)){printf("getresid: (index>=dim), %f %f %f and %d", xval, x, *weight, *index);getchar();}
}



/** INVARIANT DISTRIBUTION **/
template<size_t dim>
void inv_distri(double (&invdist)[dim], double (&prob)[dim][dim])
{
    double tempdist[dim], critdist, sumdist;
    int i, j;
    
    invdist[1] = 1.0;
    
    critdist   = 1.0;
    while(critdist > 0.00000001) {
        
        for(i=0;i<dim;i++){tempdist[i] = invdist[i];}
        
        // compute the invdist //
        for(i = 0; i<dim; i++){invdist[i] = 0.0;}
        
        for(i = 0; i<dim; i++){for(j = 0; j<dim; j++){invdist[i] += tempdist[j]*prob[j][i];}}
        
        critdist = 0.0;
        for(i=0;i<dim;i++){critdist = max(fabs(invdist[i] - tempdist[i]), critdist);}
        
//       printf("%f %f %f %f %f %f %f", invdist[0], invdist[1], invdist[2], tempdist[0], tempdist[1], tempdist[2], critdist); getchar();
        
    }
    
    // renormalize invdist //
    sumdist = 0.0;
    for(i = 0; i<dim; i++) {sumdist += invdist[i];}
    for(i=0;i<dim;i++) {invdist[i] = invdist[i]/sumdist;}
    
}







#define interPOL(x1,y1,y2) ((1.0-(x1))*(y1)+(x1)*(y2))


double GINI(double *var, double *mass, int size, double mass_tot) {
    int j,i;
    double *aera, *varcumul, *masscumul, *varcumulw, totalwealth, gini;
    
    aera = (double *) calloc((size+1), sizeof(double));
    varcumul = (double *) calloc((size+1), sizeof(double));
    masscumul = (double *) calloc((size+1), sizeof(double));
    varcumulw = (double *) calloc((size+1), sizeof(double));
    
    totalwealth = 0.0;
    masscumul[0] = 0.0;
    varcumul[0] = 0.0;
    varcumulw[0] = 0.0;
    
    for(j=0; j < size; j++) {
        totalwealth += var[j] * mass[j]/mass_tot;
    }
    
    for(j=1; j < (size+1); j++) {
        masscumul[j] = masscumul[j-1] + mass[j-1]/mass_tot;
        varcumul[j] = varcumul[j-1] + var[j-1]*mass[j-1]/mass_tot;
        varcumulw[j] = varcumul[j]/totalwealth;
    }
    
    aera[0] = 0;
    for(i=1; i < (size+1); i++) {
        aera[i] = aera[i-1] + (varcumulw[i] + varcumulw[i-1])*(masscumul[i] - masscumul[i-1]);
    }
    
    gini = 1 - aera[size];
    
    free(aera);
    free(varcumul);
    free(varcumulw);
    free(masscumul);
    
    return gini;
}





double GINI_sort(double vector[][2], int size) {
    int j,i;
    double *aera, *varcumul, *masscumul, *varcumulw, totalwealth, gini;
    
    aera = (double *) calloc((size+1), sizeof(double));
    varcumul = (double *) calloc((size+1), sizeof(double));
    masscumul = (double *) calloc((size+1), sizeof(double));
    varcumulw = (double *) calloc((size+1), sizeof(double));
    
    totalwealth = 0.0;
    masscumul[0] = 0.0;
    varcumul[0] = 0.0;
    varcumulw[0] = 0.0;
    
    for(j=0; j < size; j++) {
        totalwealth += vector[j][0] * vector[j][1];
    }
    
    for(j=1; j < (size+1); j++) {
        masscumul[j] = masscumul[j-1] + vector[j-1][1];
        varcumul[j] = varcumul[j-1] + vector[j-1][0]*vector[j-1][1];
        varcumulw[j] = varcumul[j]/totalwealth;
    }
    
    aera[0] = 0;
    for(i=1; i < (size+1); i++) {
        aera[i] = aera[i-1] + (varcumulw[i] + varcumulw[i-1])*(masscumul[i] - masscumul[i-1]);
    }
    
    gini = 1 - aera[size];
    
    free(aera);
    free(varcumul);
    free(varcumulw);
    free(masscumul);
    
    return gini;
}




// Function to get top X% of the population detening Y% of the wealth
double top_threshold_fun(double vector[][2], int size, double percent, double *dprop) {
    int iter;
    double massCUM;

    massCUM         = 0.0;
    iter            = 0;
    
    while(massCUM < percent){massCUM += vector[iter][1];iter++;}
    *dprop   = (massCUM - vector[iter-1][1]);
    //printf("%f,%f,%f,%f",percent,massCUM,dprop,inter1d(dprop,vector[iter-1][0],vector[iter][0]));getchar();
    
    return (inter1d(*dprop,vector[iter-1][0],vector[iter][0]));
}



// Function to get top X% of the population detening Y% of the wealth
double toppercent(double *var, double *mass, int size, const double percent) {
    int j,i,iter;
    double *varcumul, *varcumulw, totalwealth, massCUM, topXval, dprop;
    
    
//    printf("%f, %f, %f", percent, mass[4], var[4]); getchar();
    varcumul = (double *) calloc((size+1), sizeof(double));
    varcumulw = (double *) calloc((size+1), sizeof(double));
    
    totalwealth = 0.0;
    varcumul[0] = 0.0;
    varcumulw[0] = 0.0;
    
    for(j=0; j < size; j++) {
        totalwealth += var[j] * mass[j];
    }
    
    for(j=1; j < (size+1); j++) {
        varcumul[j] = varcumul[j-1] + var[j-1]*mass[j-1];
        varcumulw[j] = varcumul[j]/totalwealth;
    }

    massCUM = 0.0;
    iter = 0.0;
    while(massCUM < percent) {
        massCUM += mass[iter];
        iter++;
    }
    dprop = massCUM - percent;
    // Apply an 1D interpolation to get closer value to the truth.
    topXval = inter1d(dprop, varcumulw[max(0,min(size-1,iter-1))], varcumulw[max(0,min(size-1,iter))]);
//    printf("%f, %f, %f, %f, %d, %d, %f, %f, %f, %f",*topXval, dprop, inter1d(dprop, varcumulw[iter - 1], varcumulw[iter]), percent, iter, size, mass[4], var[4], varcumul[4], varcumulw[4]);getchar();
    // without inter1D:
    // topXval = varcumulw[iter - 1];
    
    free(varcumulw);
    free(varcumul);
    
    return topXval;
}






// Function to get top X% of the population detening Y% of the wealth
double toppercent_sort(double vector[][2], int size, const double percent) {
    int j,i,iter;
    double *varcumul, *varcumulw, totalwealth, massCUM, topXval, dprop;
    
    
//    printf("%f, %f, %f", percent, mass[4], var[4]); getchar();
    varcumul = (double *) calloc((size+1), sizeof(double));
    varcumulw = (double *) calloc((size+1), sizeof(double));
    
    totalwealth = 0.0;
    varcumul[0] = 0.0;
    varcumulw[0] = 0.0;
    
    for(j=0; j < size; j++) {
        totalwealth += vector[j][0] * vector[j][1];
    }
    
    for(j=1; j < (size+1); j++) {
        varcumul[j] = varcumul[j-1] + vector[j-1][0]*vector[j-1][1];
        varcumulw[j] = varcumul[j]/totalwealth;
    }

    massCUM = 0.0;
    iter = 0.0;
    while(massCUM < percent) {
        massCUM += vector[iter][1];
        iter++;
    }
    dprop = massCUM - percent;
    // Apply an 1D interpolation to get closer value to the truth.
    topXval = inter1d(dprop, varcumulw[max(0,min(size-1,iter-1))], varcumulw[max(0,min(size-1,iter))]);
//    printf("%f, %f, %f, %f, %d, %d, %f, %f, %f, %f",*topXval, dprop, inter1d(dprop, varcumulw[iter - 1], varcumulw[iter]), percent, iter, size, mass[4], var[4], varcumul[4], varcumulw[4]);getchar();
    // without inter1D:
    // topXval = varcumulw[iter - 1];
    
    free(varcumulw);
    free(varcumul);
    
    return topXval;
}








//
//  medianworth.cpp
//
//
//  Created by Alexandre GAILLARD on 01/02/2019.
//



// function to get median.
double medianworth_sort(double wealth[][2],const double fulldist, int dim)
{
    int i = 0;
    double dfrac,dfracold,prop, medianworthout;
    
    if (fulldist>0.0){
        i       =   0;
        dfrac   =   0.0;
        dfrac   +=  wealth[i][1]/fulldist;
    
        if(dfrac>=0.5){
            medianworthout=wealth[i][0];
            return(medianworthout);
        }else{
        
            while ((dfrac<0.5) && (i<(dim-1))){
                dfracold=dfrac;
                i++;
                dfrac+=wealth[i][1]/fulldist;
            }
            prop=(0.5-dfracold)/(dfrac-dfracold);
            return((1.0-prop)*wealth[(i-1)][0]+prop*wealth[i][0]);
        }
    
    }else{
//        printf("error in median");getchar();
        return(0.0);
    }
}
//







/**
 This files implement several useful functions. **/



// transform a quantile in (0,1) into normal (sigma,nu)
double norminv(double quantile, double mu, double sigma){return(mu + sigma*r8_normal_01_cdf_inverse(quantile));}

// gives the pdf of a normal distribution, for a value x
double normpdf(double mu, double sigma, double x){return(exp(-(pow((x-mu),2.0))/(2.0*pow(sigma,2.0)))/(sigma*pow((TWOPI),0.5)));}

// gives the cdf of a normal distribution, for a value x.
double normcdf(double mu, double sigma, double x){return((1.0+erf((x-mu)/(sigma*pow(2.0,0.5))))/2.0);}

// transform a quantile into a value using the pareto distribution.
double Qpareto(double quant, double xmin, double parcoef){return(xmin*pow((1.0-quant),(-1.0/parcoef)));}

// transforms normal z shock into a variable scaling wage under log-normal // pareto
double Ztrans(double z, double cdfz){
    if(cdfz < qtop){
        return(exp(z));                                        // give a value for the wage following the wage level in the log-normal
    }else{
        return(Qpareto((cdfz-qtop)/(1.0-qtop),ebar,eta));     // give a value for the wage following the wage level in the pareto
    }
}


// Gaussian Quadratures.
typedef double DP;
void gauher(double x[], double w[], int ngrid){
	const DP EPS=1.0e-14,PIM4=0.7511255444649425;
	const int MAXIT=10;
	int i,its,j,m;
	DP p1,p2,p3,pp,z,z1;

	int n = ngrid;
	m = (n+1)/2;
	for (i=0;i<m;i++) {
		if (i == 0) {z=sqrt(DP(2*n+1))-1.85575*pow(DP(2*n+1),-0.16667);
		} else if (i == 1) {z -= 1.14*pow(DP(n),0.426)/z;
		} else if (i == 2) {z=1.86*z-0.86*x[0];
		} else if (i == 3) {z=1.91*z-0.91*x[1];
		} else {z=2.0*z-x[i-2];}
		for (its=0;its<MAXIT;its++) {
			p1=PIM4;
			p2=0.0;
			for (j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=z*sqrt(2.0/(j+1))*p2-sqrt(DP(j)/(j+1))*p3;
			}
			pp=sqrt(DP(2*n))*p2;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its >= MAXIT) printf("too many iterations in gauher");
		x[i]=z;
		x[n-1-i] = -z;
		w[i]=2.0/(pp*pp);
		w[n-1-i]=w[i];
	}
}



// REDUCTION OF A TRANSITION MATRIX
// This function transform a transition matrix of P[M][M'] into a lower scale matrix. For each M, we save the number of
// index that are non-zero, inx_pos[M][max_P] save the 'non-zero' index values.
void reduction_matrix(){
    
    int n,n2,n3,i,ii;
    int size = length_y;
    double trunc=1e-06, sum_P = 0;
    
    for(i=0;i<size;i++){
        sum_P = 0;
        for(ii=0;ii<size;ii++){sum_P      += P_Y[i][ii];}
        for(ii=0;ii<size;ii++){P_Y[i][ii]  = P_Y[i][ii]/sum_P;}
        
        n3 = 0;
        for(ii=0;ii<size;ii++){
            if(P_Y[i][ii] > trunc){
                P_Y_reduc[i][n3]  = ii;
                n3                = n3 + 1;
            }else{
                P_Y[i][ii]  = 0.0;
            }
        }
        
        sum_P = 0;
        for(ii=0;ii<size;ii++){sum_P     += P_Y[i][ii];}
        for(ii=0;ii<size;ii++){P_Y[i][ii] = P_Y[i][ii]/sum_P;}
        
        max_P_Y[i] = n3;   //gives the index for i.
    }
}




//  This function gives the weight of each point in the grid of EEGRID, using a linear weighting rule.
template<size_t Nodes2>
void intWeightsN(double mu, double sigma, double zzgrid[Nodes2], double eegrid[Nodes2], double (&wx)[Nodes2]){
    
    int ipts = 100000;
    double intgridZ[ipts] = {0};
    double extrap = 0.1;        // extrapolate because can be off-the-grid.
    
    int ii,n, zz,y;
    double pdfval,weight,xlo,xhi,eval;
    
    xlo = (1.0+extrap)*zzgrid[0] - extrap*zzgrid[Nodes2-1];
    xhi = (1.0+extrap)*zzgrid[Nodes2-1] - extrap*zzgrid[0];
    
    for(ii=0; ii<Nodes2; ii++){wx[ii] = 0.0;}
    
    for(ii=0; ii<ipts; ii++){
      intgridZ[ii]  = linspace(xlo,xhi,ipts,ii);
      eval          = Ztrans(intgridZ[ii],normcdf(mu,sigma,intgridZ[ii])); // cdf.
      
      // locate eval on the grid eegrid.
      zz = 0;
      while(eval > eegrid[zz]){zz = zz + 1;}
      zz = zz-1;
      n         = max(0,min(Nodes2-2,zz));
      weight    = (eval-eegrid[zz])/(eegrid[zz+1]-eegrid[zz]);
      
      pdfval    = normpdf(mu,sigma,intgridZ[ii]); // probability
      wx[n]     += (1.0-weight)*pdfval;
      wx[n+1]   += weight*pdfval;
      
//      printf("%d | %f %f, %f %f %d %f | %f | %f %f", ii, xlo, xhi, intgridZ[ii], eval, n, weight, pdfval, wx[n], wx[n+1]);getchar();
    }
    
    double sum_W = 0;
    for(y=0;y<length_y;y++){sum_W += wx[y];}
    for(y=0;y<length_y;y++){wx[y] = wx[y]/sum_W;}
    
}

