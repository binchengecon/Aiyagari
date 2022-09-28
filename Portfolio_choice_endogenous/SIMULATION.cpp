//
//  simulation.cpp
//  Created by Alexandre GAILLARD on 30/01/2019.
//

/* ALGORITHM: in order to simulate the model, I proceed as follows:
    - I use the non-stochastic method developed in Young (2010) -- (with improvements)
    - This part is completely invariante from aggregate dynamics.       */


void SIMULATION(double *DIST,double *save, double *cons, double *NU, double *X, double *N, int transi_dyn){

    
// initialize the distribution + indexes
int     i,yt,y,t,ii,z,ipsi,ynext,tnext;
double  psi_next,s_next,weight,distverif;


double  *DIST_OLD;  DIST_OLD = (double *) calloc((length_x), sizeof(double));
double  *DIST_TMP;  DIST_TMP = (double *) calloc((length_x), sizeof(double));
FILE    *dfile;
   
   
//makezero(dist,length_x);
//dist[0] = 1.0;


// START LOOP -- SIMULATION
double CRIT_DIST = 1.0; int iter_simul = 0; int stop_transi = 0;
while(CRIT_DIST > EPS_DIST && (iter_simul < MAX_ITER_SIM) && (stop_transi == 0)){


    // copy distin in distold //
    bascule(DIST,DIST_OLD,length_x);
    makezero(DIST,length_x);
    makezero(DIST_TMP,length_x);
    
    CRIT_DIST       = 0.0;  distverif        = 0.0;
    
    // compute the distribution next
    #if OPT_OMP == 1
        #pragma omp parallel for reduction(+:DIST_TMP[:length_x]) private(yt,y,t,i,z,psi_next,ipsi,weight)
    #endif
    for(yt = 0; yt < length_yt; yt++){      
    
        y = (int)(floor((yt/(length_t))));       // find index y
        t = (int)(yt-length_t*y);                // find index t
            
        for(i = 0; i < length_i; i++){     // this is the end-of-period decision of saving' (how much to save)
            if(DIST_OLD[x(i,y,t)] > 0){    // check if distribution exists.
                for(z = 0; z < length_z; z++){
                
                    psi_next    = (PSI[i]+W*Y[y]) + NU[x(i,y,t)]*(PSI[i]+W*Y[y])*(Z[z]+premium) + (1.0-NU[x(i,y,t)])*(PSI[i]+W*Y[y])*(0.04);
                    
                    ipsi     = max(0,min((int)(floor(invgridspace(psi_next,PSI_min,PSI_max,length_i))),length_i-2)); // (also defined on igrid)
                    weight   = min(1.0,(psi_next - PSI[ipsi])/(PSI[ipsi+1] - PSI[ipsi]));
                    if(ipsi < 0 | weight < 0.0 | weight > 1.0){printf("mistake 1: %d %f psi_next = %f", ipsi, weight, psi_next);getchar();}
                    
                    // distribute to the temp distribution.
                    DIST_TMP[x(ipsi,y,t)]      += (1.0-weight)*TZ[z]*DIST_OLD[x(i,y,t)];
                    DIST_TMP[x((ipsi+1),y,t)]  += weight*TZ[z]*DIST_OLD[x(i,y,t)];
                    
                } // end of e
            } // end of DIST_OLD
        } // end of i
    } // end of yt
    
    
    #if OPT_OMP == 1
        #pragma omp parallel for reduction(+:DIST[:length_x]) private(yt,y,t,i,ynext,tnext,s_next,ipsi,weight)
    #endif
    for(yt = 0; yt < length_yt; yt++){
    
        y = (int)(floor((yt/(length_t))));       // find index y
        t = (int)(yt-length_t*y);                // find index t
            
        for(i = 0; i < length_i; i++){     // this is the end-of-period decision of saving' (how much to save)
            for(ynext = 0; ynext < length_y; ynext++){
                for(tnext = 0; tnext < length_t; tnext++){
                
                    // if don't die.
                    s_next  = (save[x(i,y,t)])*(1.0-TAX_wealth);
                    ipsi     = max(0,min((int)(floor(invgridspace(s_next,A_min,A_max,length_i))),length_i-2)); // (also defined on igrid)
                    weight   = min(1.0,(s_next - A[ipsi])/(A[ipsi+1] - A[ipsi]));
                    if(ipsi < 0 | weight < 0.0 | weight > 1.0){printf("mistake 1: %d %f psi_next = %f", ipsi, weight, psi_next);getchar();}
                                        
                    DIST[x(ipsi,ynext,tnext)]       += (1.0-weight)*P_Y[y][ynext]*P_THETA[t][tnext]*(1.0-P_die)*DIST_TMP[x(i,y,t)];
                    DIST[x((ipsi+1),ynext,tnext)]   += weight*P_Y[y][ynext]*P_THETA[t][tnext]*(1.0-P_die)*DIST_TMP[x(i,y,t)];
                    
                    // if die.
                    s_next  = (save[x(i,y,t)])*(1.0-TAX_bequest);
                    ipsi     = max(0,min((int)(floor(invgridspace(s_next,A_min,A_max,length_i))),length_i-2)); // (also defined on igrid)
                    weight   = min(1.0,(s_next - A[ipsi])/(A[ipsi+1] - A[ipsi]));
                    if(ipsi < 0 | weight < 0.0 | weight > 1.0){printf("mistake 1: %d %f psi_next = %f", ipsi, weight, psi_next);getchar();}
                                        
                    DIST[x(ipsi,ynext,0)]       += (1.0-weight)*P_Y[y][ynext]*P_THETA[t][tnext]*(P_die)*DIST_TMP[x(i,y,t)];
                    DIST[x((ipsi+1),ynext,0)]   += weight*P_Y[y][ynext]*P_THETA[t][tnext]*(P_die)*DIST_TMP[x(i,y,t)];
                    
                }  // end tnext
            } // end ynext
        } // end of i
    } // end of yt
    
    
    iter_simul++;
    
    // convergence criterion
    if(iter_simul < MAX_ITER_SIM && transi_dyn == 0){
    
        #if OMP == 1
        #pragma omp parallel for reduction(+:verif) reduction(max:CRIT_DIST) private(ii)
        #endif
        for(ii = 0; ii < length_x; ii++){      // productivity wage
            CRIT_DIST    = max(CRIT_DIST,fabs(DIST[ii]-DIST_OLD[ii]));
            distverif   += DIST[ii];
        } // end y
    }
    
    if(transi_dyn == 1){stop_transi = 1;}
    
//   printf("CNVG SIMUL = %f | distverif = %f \n", CRIT_DIST, distverif);

}   // end while for simulation





// check for non-stationarity + ergodicity.
if(CRIT_DIST > EPS_DIST && transi_dyn == 0){printf("\n \n possibly no stationarity :: cnvg = %f %f %d \n \n",CRIT_DIST,distverif,iter_simul);}



// ---- STATISTICS: compute the relevant statistics.

// put to zero and initialize.
prop_investor_model = 0.0;
Q_tot = 0.0; L_tot = 0; K_tot = 0.0;

// several distribution in function of wealth.
// WEIGHT gives the distribution of weight for a wealth level "i".
for(i = 0; i < length_i; i++){
    WEIGHT[i]                       = 0.0;
    DIST_riskyasset[i]              = 0.0;
    DIST_fracinvestor[i]            = 0.0;
    DIST_productivity[i]            = 0.0;
    DIST_returns_conditional[i]     = 0.0;
    DIST_returns_unconditional[i]   = 0.0;
    DIST_earnings[i]                = 0.0;
}

double Q_E = 0.0, Q_c = 0.0, Q_debt = 0.0;
for(i = 0; i < length_i; i++){
    for(yt = 0; yt < length_yt; yt++){      // next period riskfree asset
        
        y = (int)(floor((yt/(length_t))));       // find index y
        t = (int)(yt-length_t*y);                // find index t
        
        WEIGHT[i]            += DIST[x(i,y,t)];                 // weight for wealth distribution.
        DIST_riskyasset[i]   += DIST[x(i,y,t)]*NU[x(i,y,t)];    // distribution of portfolio choice.
        DIST_productivity[i] += DIST[x(i,y,t)]*THETA[t];        // distribution productivity
        DIST_earnings[i]     += DIST[x(i,y,t)]*W*Y[y];    // distribution earnings
        DIST_returns_unconditional[i]  += DIST[x(i,y,t)]*(NU[x(i,y,t)]*(PSI[i]+W*Y[y])*(Z[z]+premium) + (1.0-NU[x(i,y,t)])*(PSI[i]+W*Y[y])*(0.04))/(PSI[i]+W*Y[y]);        // distribution returns all wealth
//        if(shareNU > 0.0){DIST_returns_conditional[i]  += DIST[x(i,y,t)]*(THETA[t]*RQ*profit(K_val) - d_borr*R - (deltapar)*K_val)/(shareNU*A[i]);}
        // distribution returns private equity
        
        // distribution investors
        if(NU[x(i,y,t)] > 0.0001){DIST_fracinvestor[i] += DIST[x(i,y,t)];}
        if(NU[x(i,y,t)] > 0.0001){prop_investor_model  += DIST[x(i,y,t)];}
    
        // distribution disaggregated.
//        DIST_earnings_x[x(i,y,t)][0]   = (1.0-fun_alloc_time(shareNU))*Y[y];
//        DIST_earnings_x[x(i,y,t)][1]   = DIST[x(i,y,t)];
        DIST_riskyasset_x[x(i,y,t)][0] = NU[x(i,y,t)]*PSI[i];
        DIST_riskyasset_x[x(i,y,t)][1] = DIST[x(i,y,t)];
        
        // aggregate quantities
        L_tot   += DIST[x(i,y,t)]*Y[y];
        K_tot   += DIST[x(i,y,t)]*(PSI[i]);
        for(z=0; z<length_z; z++){Q_tot  += TZ[z]*DIST[x(i,y,t)]*(NU[x(i,y,t)]*Z[z]*PSI[i] + (1.0-NU[x(i,y,t)])*PSI[i]*TFP);}
        
    }
    

    // RENORMALIZE BY WEIGHT
    DIST_riskyasset[i]             = DIST_riskyasset[i]/WEIGHT[i];
    DIST_productivity[i]           = DIST_productivity[i]/WEIGHT[i];
    DIST_earnings[i]               = DIST_earnings[i]/WEIGHT[i];    // distribution of persistent earnings.
    DIST_returns_unconditional[i]  = DIST_returns_unconditional[i]/WEIGHT[i];
    DIST_returns_conditional[i]    = DIST_returns_conditional[i]/WEIGHT[i];
    
    
}

//Q_tot = Q_E + Q_c - Q_debt;
GDP = production(Q_tot,L_tot);

// qsort
qsort(DIST_earnings_x,(length_x),2*sizeof(double),comparefun2);
qsort(DIST_riskyasset_x,(length_x),(2)*sizeof(double),comparefun2);

// compute Gini coefficient;
GINI_wealth      =   GINI(PSI, WEIGHT, length_i, 1.0);
GINI_riskyasset  =   GINI_sort(DIST_riskyasset_x,(length_x));
GINI_earnings    =   GINI_sort(DIST_earnings_x,(length_x));

// compute the share in different quantiles of wealth.
top50       =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.50);
top60       =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.60);
top70       =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.70);
top80       =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.80);
top90       =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.90);
top95       =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.95);
top99       =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.99);
top99_9     =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.999);
top99_99    =   1.0 - toppercent(PSI, WEIGHT, length_i, 0.9999);

// compute the share in different quantiles of risky assets.
top50_riskasset        =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.50);
top60_riskasset        =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.60);
top70_riskasset        =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.70);
top80_riskasset        =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.80);
top90_riskasset        =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.90);
top95_riskasset        =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.95);
top99_riskasset        =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.99);
top99_9_riskasset      =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.999);
top99_99_riskasset     =   1.0 - toppercent_sort(DIST_riskyasset_x,(length_x), 0.9999);


/* SET THE MOMENT OF THE MODEL */
set_genmoments[0] = K_tot/GDP;
set_genmoments[1] = Q_tot/GDP;
set_genmoments[2] = prop_investor_model;
set_genmoments[3] = top99;
set_genmoments[4] = top95_riskasset;



/*******************************************************/
/// --- PRINT VARIOUS RESULTS

//if(show_stats == 1){
FILE *logfile;logfile = fopen(FILE_log, "a");

fprintf(logfile,"------ Estimated Model parameters ------\n");
printf("------ Estimated Model parameters ------\n");
for(int p=0; p<nb_para; p++){fprintf(logfile,"%s = %f\n",set_params_name[p],set_params[p]);}
for(int p=0; p<nb_para; p++){printf("%s = %f\n",set_params_name[p],set_params[p]);}
fprintf(logfile,"\n");

fprintf(logfile,"------ Estimated Moment [target] ------\n");
printf("------ Estimated Moment [target] ------\n");
for(int m=0; m<nb_moments; m++){fprintf(logfile,"%s \t obs = %f, gen = %f \n",set_moments_name[m],set_obsmoments[m],set_genmoments[m]);}
for(int m=0; m<nb_moments; m++){printf("%s \t obs = %f, gen = %f \n",set_moments_name[m],set_obsmoments[m],set_genmoments[m]);}
fprintf(logfile,"\n");

fprintf(logfile,"------ STATISTICS ------\n");
printf("------ GINI COEFFICIENT------\n");
fprintf(logfile,"GINIearnings = %f   |  GINIwealth = %f   |  GINIriskyassets = %f \n", GINI_earnings, GINI_wealth, GINI_riskyasset);
fprintf(logfile,"------ TOP SHARES ------\n");
fprintf(logfile,"wealth, [50, 60, 70, 80, 90, 95, 99, 99.9, 99.99] = [%f, %f, %f, %f, %f, %f, %f, %f, %f] \n",top50,top60,top70,top80,top90,top95,top99,top99_9,top99_99);
fprintf(logfile,"risky assets, [50, 60, 70, 80, 90, 95, 99, 99.9, 99.99] = [%f, %f, %f, %f, %f, %f, %f, %f, %f] \n",top50_riskasset,top60_riskasset,top70_riskasset,top80_riskasset,top90_riskasset,top95_riskasset,top99_riskasset,top99_9_riskasset,top99_99_riskasset);
printf("GINIearnings = %f   |  GINIwealth = %f   |  GINIriskyassets = %f \n", GINI_earnings, GINI_wealth, GINI_riskyasset);
printf("------ TOP SHARES ------\n");
printf("wealth, [50, 60, 70, 80, 90, 95, 99, 99.9, 99.99] = [%f, %f, %f, %f, %f, %f, %f, %f, %f] \n",top50,top60,top70,top80,top90,top95,top99,top99_9,top99_99);
printf("risky assets, [50, 60, 70, 80, 90, 95, 99, 99.9, 99.99] = [%f, %f, %f, %f, %f, %f, %f, %f, %f] \n",top50_riskasset,top60_riskasset,top70_riskasset,top80_riskasset,top90_riskasset,top95_riskasset,top99_riskasset,top99_9_riskasset,top99_99_riskasset);

fprintf(logfile,"------ EQUILIBRIUM VARIABLES ------\n");
printf("------ Endogenous Prices ------\n");
fprintf(logfile,"Prices: [R,W] = [%20.15f, %20.15f] \n",R,W);
fprintf(logfile,"Q = %20.15f, L = %20.15f, K = %20.15f, Y = %20.15f, K/Y = %20.15f \n",Q_tot,L_tot,K_tot,GDP,K_tot/GDP);
printf("Prices: [R,W] = [%20.15f, %20.15f] \n",R,W);
printf("Q = %20.15f, L = %20.15f, K = %20.15f, Y = %20.15f, K/Y = %20.15f \n",Q_tot,L_tot,K_tot,GDP,K_tot/GDP);



#if PRINT_SS == 1 | PRINT_SS == 2

    // POLICY FUNCTIONS //
    dfile=fopen(FILE_dist, "w"); setbuf (dfile, NULL );
    for(i = 0; i < length_i; i++){
        fprintf(dfile,"%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",A[i],WEIGHT[i],DIST_riskyasset[i],DIST_fracinvestor[i],DIST_productivity[i],DIST_earnings[i],DIST_returns_conditional[i],DIST_returns_unconditional[i]);
    } // end i
    fclose(dfile);
    
#endif



// FREE POINTERS
free(DIST_OLD);
free(DIST_TMP);


}

