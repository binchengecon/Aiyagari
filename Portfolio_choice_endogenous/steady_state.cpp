/** -------------------------------------------------
    Alexandre GAILLARD - PHILIPP WANGNER  (GIF PAPER, 2020)
    
    -------
    STEADY_STATE ::
    This file compute the long-run steady state of the model.
        - OUTER_LOOP_aggQ   :: compute the resulting aggregates given a Q value.
        - steady_state      :: compute the (long-run) steady-state of the model.
    -------
**/
    

/** -------------------- GDP, PRICES -------------------- **/
void price_function(const double AGG_Q, const double AGG_L){
    GDP = production(AGG_Q,AGG_L);
    RQ  = MPQ(AGG_Q,AGG_L); 
    R   = MPQ(AGG_Q,AGG_L); //risk_free;
    W   = MPL(AGG_Q,AGG_L);
}


// initialize struct
struct structQ_OUTER{double *save_struct,*cons_struct,*X_struct,*N_struct,*NU_struct,*VF_struct,*dVF_struct,*DIST_struct,Qnew_struct;int iter_Q_struct;};


// OUTER-LOOP :: Function to compute aggregate capital + labor
double OUTER_LOOP_aggQ(const double AGG_Q, void *params){

    // initialize timer.
    timeval t1, t2, t3, t4;
    
    struct structQ_OUTER *QPARAMS = (struct structQ_OUTER *) params;
    int    iter_Q    = QPARAMS->iter_Q_struct;
    // policy function
    double *save     = QPARAMS->save_struct;
    double *cons     = QPARAMS->cons_struct;
    double *NU       = QPARAMS->NU_struct;
    double *X        = QPARAMS->X_struct;
    double *N        = QPARAMS->N_struct;
    // value function
    double *VF       = QPARAMS->VF_struct;
    double *dVF      = QPARAMS->dVF_struct;
    // distribution
    double *DIST     = QPARAMS->DIST_struct;


    // compute prices
    price_function(AGG_Q,L_tot);
    printf("- OUTER LOOP PRICES: Qold = %f, R = %f, RQ = %f,  W = %f, L_tot = %f \n", AGG_Q, R, RQ, W, L_tot);


    // policy functions
    gettimeofday(&t1, NULL);
        POLICY_EGM(VF,dVF,save,cons,NU,X,N,0);
    gettimeofday(&t2, NULL);
    printf("-- POLICY DONE in %fs \n", timer_fun(t1,t2));
    printf("-- \n");


    // simulation
    gettimeofday(&t3, NULL);
//        SIMULATION(DIST,save,cons,NU,X,N,0);
    gettimeofday(&t4, NULL);
    printf("-- \n");
    printf("-- SIMULATION DONE in %fs \n", timer_fun(t3,t4));

    
    // Q_tot is computed within SIMULATION/
    QPARAMS->iter_Q_struct  =   iter_Q+1;
    QPARAMS->Qnew_struct    =   Q_tot;

    // report GDP + statistics.
    GDP = production(Q_tot,L_tot);
//    printf("-- STATS: Q = %f, Qnext = %f, dist = %f, K = %f, Y = %f, K/Y = %f, Gini = [W=%f, Risky=%f, E=%f], top = [%f %f %f], R = %f \n", AGG_Q, Q_tot, AGG_Q - Q_tot, K_tot, GDP, K_tot/GDP, GINI_wealth, GINI_riskyasset, GINI_earnings, top90, top99, top99_9, prop_investor_model, R);
    
    return (AGG_Q - Q_tot);
}





// STEADY-STATE
double steady_state(double *set_params){


    // change the set of parameters + processes.
    set_params_eq(set_params);
    set_processes();


    /** ------------------ INITIALIZING ------------------ **/
    
    // initialize integer
    int	a, i,y,t,k,z, iter;
    
    // initialize timer
    timeval t1, t2, t3, t4;
    double elapsedTime;
    
    // pointer of the model
    double *VF, *dVF, *save, *cons, *X, *NU, *N, *DIST;
    VF      = (double *) calloc((length_x), sizeof(double));
    dVF     = (double *) calloc((length_x), sizeof(double));
    NU      = (double *) calloc((length_x), sizeof(double));
    DIST    = (double *) calloc((length_x), sizeof(double));
    save    = (double *) calloc((length_x), sizeof(double));
    cons    = (double *) calloc((length_x), sizeof(double));
    X       = (double *) calloc((length_x), sizeof(double));
    N       = (double *) calloc((length_x), sizeof(double));


    // guess for prices.
    R       = R_guess;
    W       = W_guess;
    L_tot   = L_guess;

    // guess for value function.
    for(y = 0; y < length_y; y++){
        for(t = 0; t < length_t; t++){
            for(i = 0; i < length_i; i++){
                VF[x(i,y,t)] = U(PSI[i]*(1.0+R) + W*Y[y])/(1.0-beta);
                derivative(i,y,t,VF,dVF);
            }
        }
    }
    
    // initial distribution for steady-state (for transition dynamic, gonna be different).
    DIST[0] = 1.0;
    


    /** ------------------ EQUILIBRIUM ------------------
        - The routine use a bisection to find the equilibrium Q.
    **/
    
    printf("STARTING EQUILIBRIUM...\n");


    double Q_new = (EQ_QMAX + EQ_QMIN)/2.0;
    double CRIT_equilibrium = 1.0;
    while(CRIT_equilibrium > EPS_equilibrium){

        struct structQ_OUTER paramsQ;
        paramsQ.save_struct         = save;
        paramsQ.cons_struct         = cons;
        paramsQ.N_struct            = N;
        paramsQ.X_struct            = X;
        paramsQ.NU_struct           = NU;
        paramsQ.VF_struct           = VF;
        paramsQ.dVF_struct          = dVF;
        paramsQ.DIST_struct         = DIST;
        paramsQ.iter_Q_struct       = 0;
        
        OUTER_LOOP_aggQ(Q_new,&paramsQ);

        CRIT_equilibrium = fabs(Q_new - Q_tot);
        
        Q_new = 0.01*Q_tot + 0.99*Q_new;
        
        printf("\n \n CRIT_EQUILIBIRUM = %f \n \n", CRIT_equilibrium);
    }




    /** ------------------ COMPUTE THE MOMENTS OF THE MODEL ------------------
        - indirect inference on the moments using the simulated method of moments.
    **/

    // Step 1:: characterize the covariance-matrix of moments.
    double covmat[nb_moments][nb_moments], SMM_dist;
    for(int m = 0; m < nb_moments; m++){covmat[m][m] = 1.0;}    // put weight of 1 to the diagonal matrix

    // Step 2:: compute the distance induced by the observed and the simulated moments.
    SMM_dist = SMM(set_obsmoments,set_genmoments,covmat);
    //if(stop_GMM == 1){SMM_dist = 10000;}        // security.

    printf("- SMM-RESULT: Value = %f, TIME = %fs \n", SMM_dist, timer_fun(t3,t4));

    FILE *logfile;logfile = fopen(FILE_log, "a");
        fprintf(logfile,"- SMM-RESULT: Value = %f, TIME = %fs \n\n", SMM_dist, timer_fun(t3,t4));
    fclose(logfile);


    // free up the pointers.
    free(VF);
    free(dVF);
    free(save);
    free(cons);
    free(X);
    free(NU);
    free(N);
    free(DIST);

    // return the value of the SMM distance
    return SMM_dist;

}


