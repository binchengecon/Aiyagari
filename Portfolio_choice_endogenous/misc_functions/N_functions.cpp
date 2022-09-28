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


/**
 This files implement several useful functions. **/

// useful numbers:
const double PI       = 3.141592653589793238462643383279502884197;
const double PIO2     = 1.57079632679489661923132169163975144209858;
const double TWOPI    = 6.283185307179586476925286766559005768394;
const double SQRT2    = 1.41421356237309504880168872420969807856967;
const double EULER    = 0.5772156649015328606065120900824024310422;


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
        return(exp(z));                                         // give a value for the wage following the wage level in the log-normal
    }else{
        return(Qpareto((cdfz-qtop))/((1.0-qtop),ebar,eta));     // give a value for the wage following the wage level in the pareto
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
		if (its >= MAXIT) nrerror("too many iterations in gauher");
		x[i]=z;
		x[n-1-i] = -z;
		w[i]=2.0/(pp*pp);
		w[n-1-i]=w[i];
	}
}



// REDUCTION OF A TRANSITION MATRIX
// This function transform a transition matrix of P[M][M'] into a lower scale matrix. For each M, we save the number of
// index that are non-zero, inx_pos[M][max_P] save the 'non-zero' index values.
void reduction_matrix(double P_TRANS[][],int size,int inx_pos[][],int max_P){
    
    int n,n2,n3;
    double trunc=1e-06, sum_P = 0;
    
    for(i=0;i<size;i++){
        sum_P = 0;
        for(ii=0;ii<size;ii++){sum_P += P_TRANS[i][ii];}
        for(ii=0;ii<size;ii++){P_TRANS[i][ii] += P_TRANS[i][ii]/sum_P;}
        
        n3 = 0;
        for(ii=0;ii<size;ii++){
            if(P_TRANS[i][ii] > trunc){
                inx_pos[i][n3]  = ii;
                n3              = n3 + 1;
            }else{
                P_TRANS[i][ii]  = 0.0;
            }
        }
        
        sum_P = 0;
        for(ii=0;ii<size;ii++){sum_P          += P_TRANS[i][ii];}
        for(ii=0;ii<size;ii++){P_TRANS[i][ii] += P_TRANS[i][ii]/sum_P;}
        
        max_P[i] = n3+1;   //gives the index for i.
    }
}




//  This function gives the weight of each point in the grid of EEGRID, using a linear weighting rule.
double intWeightsN(double mu, double sigma, double zzgrid, double eegrid, int ngrid, double wx[]){
    
    int ipts = 100000;
    int intgridZ[ipts] = {0};
    double extrap = 0.1;        // extrapolate because can be off-the-grid.
    
    int ii,n, zz;
    double pdfval,weight,xlo,xhi,eval;
    
    xlo = (1.0+extrap)*zzgrid[0] - extrap*zzgrid[ngrid-1];
    xhi = (1.0+extrap)*zzgrid[ngrid-1] - extrap*zzgrid[0];
    
    for(ii=0; ii<ipts; ii++){
      intgridZ[ii]  = linspace(xlo,xhi,ipts,ii);
      eval          = Ztrans(intgridZ[ii],normcdf(mu,sigma,intgridZ[ii])); // cdf.
      
      // locate eval on the grid eegrid.
      zz = 0;
      while(eval > eegrid[zz]){zz = zz + 1;}
      zz = zz-1;
      n         = max(1,min(ngrid-2,zz));
      weight    = (eval-eegrid[zz])/(eegrid[zz+1]-eegrid[zz});
      
      pdfval    = normpdf(mu,sigma,intgridZ[ii]); // probability
      wx[n]     = wx[n] + (1.0-weight)*pdfval;
      wx[n+1]   = wx[n+1] + weight*pdfval;
    }
    
    double sum_W = 0;
    for(y=0;y<length_y;y++){sum_W += wx[y];}
    for(y=0;y<length_y;y++){wx[y] = wx[y]/sum_W;}
    
}

