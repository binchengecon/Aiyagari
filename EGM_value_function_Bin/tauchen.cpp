
// Original Contributor: SUMUDU KANKANAMGE (2017)
// Current Contributor: Bin (2022)

// Purpose: 
// Approximate an AR(1) process with a N-states (node) first order process

// Detail:
// Suppose we have AR(1) as y_t = (1-rho)*mu_y + rho*y_{t-1}+ epsilon_t, where epsilon_t is N(0,sigma^2) 
// We use n nodes to approxiamate this continuously valued markov chain
// The space Y={y_0,y_2,..,y_{n-1}}
// The transition matrix is P={p_{i,j}}_{0<= i,j <n}
// 1. sigma_y = sigma/sqrt(1-rho^2)
// 2. y_0 = mu_y - m*sigma_y, y_{n-1}=mu_y+m*sigma_y, delta_y = (y_{n-1}-y_0)/(n-1)
// 3.                               0=j,                    (1-rho)*mu_y + rho*y_i+epsilon <  y_j+delta_y/2
//              p_{i,j} = {         0<j<n-1, y_j-delta_y/2 < (1-rho)*mu_y + rho*y_i+epsilon <  y_j+delta_y/2 }
//                                  n-1=j,                  y_j-delta_y/2 < (1-rho)*mu_y + rho*y_i+epsilon



//--------Functions definition:

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




template<size_t Nodes>
void tauchenfun(double rho, double m, double amu, double sigma, double grid[Nodes], double Ptransi[Nodes][Nodes])

{

	
	
	int i,j,k;
	
	
	
	//First lets compute unconditiontal variance of yt
	
	double varyt= (sigma*sigma)/(1-rho*rho);
	
	//Compute stddev of yt
	
	double stdyt= sqrt(varyt);
	// std::cout << "value of std_e="<<stdyt <<"\n";
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

