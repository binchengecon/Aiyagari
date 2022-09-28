/** -------------------------------------------------
    Alexandre GAILLARD - PHILIPP WANGNER  (GIF PAPER, 2020)
    
    -------
    CALIBRATION ::
    This file specify the parameters of the model that are used to match the data moments of interest.
        - specify_params    :: specify the parameters of the model, changing endogeneously
        - specify_moments   :: specify the moments of the model, fixed
        - derivative        :: to compute the numerical derivative of the value function.
        - set_process       :: set the process of the model: labor income, investment etc.
    -------
**/


/** -------------------- ENDOGENOUS PARAMETER --------------------
    Here change
        -   specify_params + set_params_eq :: the endogenous parameter of the model.
        -   specify_moments :: the moments of interest
    
**/
void specify_params() {
    // set the parameter mean.
    set_params[0] = 0.93;   strncpy(set_params_name[0], "Beta", 10);                     // Beta - discount rate
    set_params[1] = 1.0;   strncpy(set_params_name[0], "DRS", 10);                       // DRS  - decreasing return to scale scale
    set_params[2] = 0.05;  strncpy(set_params_name[0], "scale of risk aversion", 10);   // gammapar0
    set_params[3] = 1.0;   strncpy(set_params_name[0], "elasticity parameter", 10);     // gammapar1
    
    // specify the bound for SMM
    for(int p=0;p<nb_para;p++){set_params_min[p] = set_params[p]*(0.9); set_params_max[p] = set_params[p]*1.1;}
}

void set_params_eq(double *set_params) {
    beta            = set_params[0];
    mupar           = set_params[1];
//    gammapar0       = set_params[2];
//    gammapar1       = set_params[3];
}

void specify_moments() {
    set_obsmoments[0] = 3.20;   strncpy(set_moments_name[0], "K/Y ratio/10", 100);
    set_obsmoments[1] = 0.20;   strncpy(set_moments_name[1], "Calibrated to match decreasing pattern of returns to PE", 100);
    set_obsmoments[2] = 0.13;   strncpy(set_moments_name[2], "Mean of risk-taking", 100);
    set_obsmoments[3] = 0.30;   strncpy(set_moments_name[3], "Elasticity of returns to wealth", 100);
}




/** -------------------- COMPUTE DERIVATIVES -------------------- **/
void derivative(int i, int y, int t, double *value, double *dvalue){
    if(i>=2){dvalue[x((i-1),y,t)] = deriv(PSI[i-2],PSI[i-1],PSI[i],value[x((i-2),y,t)],value[x((i-1),y,t)],value[x(i,y,t)]);}
    if(i==(length_i-1)){dvalue[x(i,y,t)] = (value[x(i,y,t)]-value[x((i-1),y,t)])/(PSI[i]-PSI[(i-1)]);}
    if(i == 1){dvalue[x((i-1),y,t)] = (value[x(1,y,t)]-value[x(0,y,t)])/(PSI[1]-PSI[0]);}
}



/** -------------------- SET THE PROCESS OF THE MODEL --------------------
    We specify the different processes;
        -   Earning process,
        -
**/
void set_processes() {


    // initialize various quantities.
    int	j,a,i,y,t,k,z,nu,iter,yy,e,jj,yinx;
    double sum_PY = 0.0, sum_L = 0.0, ynext_tmp, ynext_val, yweight;
    
    

    /** ------------------ EARNING PROCESS ------------------ **/
    
    // Persistent process for earnings.
    #if OPT_yprocess == 1
        if(length_y == 1){
            Y[0] = 1.0; P_Y[0][0] = 1.0;
        }else{

            // initialize
            double Y_tmp[length_y] = {0}, P_Y_tmp[length_y] = {0}, ExGQ[length_y] = {0}, EwGQ[length_y] = {0};
            
            //  unconditional grid and weights
            double std_y2 = std_y/sqrt(1.0-pow(rho_y,2.0));         // uncond. st. dev. persistent earnings shock
            double mu_y   = -(pow(std_y2,2.0))/2.0;                 // autocorrelation persistent earnings shock
            double zbar   = norminv(qtop,mu_y,std_y2);              // maximum value of the log-normal distribution.
            ebar   = exp(zbar);                              // rescale for the state value (min value for the pareto law)
            
            for(y = 0; y < length_y; y++){
                Y_tmp[y] = norminv(Ugrid[y],mu_y,std_y2);           // log-normal.
                Y[y]     = Ztrans(Y_tmp[y],Ugrid[y]);               // split between log-normal, or pareto.
//                printf("%f %f %f\t",Y_tmp[y],Y[y],Ugrid[y]); getchar();
            }
            
            // compute the weight for each level of "y" (should sum to one).
            intWeightsN(mu_y,std_y2,Y_tmp,Y,I_Y);
            for(y = 0; y<length_y; y++){printf("%f, %f\n",Y[y],I_Y[y]);}printf("\n");
            

            //  conditional grid and weights
            gauher(ExGQ,EwGQ,length_y2);   // weights for the gaussian hermite quadrature, notice that we consider only 10 states for persistent shocks.
            sum_PY = 0;
            for(j = 0; j < length_y2; j++){sum_PY += EwGQ[j];}
            for(j = 0; j < length_y2; j++){EwGQ[j] = EwGQ[j]/sum_PY;}

            // generate the transition matrix between states.
            for(j = 0; j < length_y; j++){for(jj = 0; jj < length_y2; jj++){P_Y[j][jj] = 0.0;}}
            for(j = 0; j < length_y; j++){
                for(jj = 0; jj < length_y2; jj++){
                    ynext_tmp = rho_y*Y_tmp[j] + (1-rho_y)*mu_y + std_y*ExGQ[jj]*sqrt(2.0);   // compute the next value given Zgrid[j]
                    ynext_val = Ztrans(ynext_tmp,normcdf(mu_y,std_y2,ynext_tmp));             // give the value of Zval, provided its cdf and value notice that for the pareto, only the cdf matters.

                    // locate ynext_val on the grid eegrid.
                    yy = 0; while(ynext_val > Y[yy]){yy++;}
                    yy      = yy-1;
                    yinx    = max(0,min(length_y-2,yy));
                    yweight = (ynext_val-Y[yinx])/(Y[yinx+1]-Y[yinx]);
                    P_Y[j][yinx]    += EwGQ[jj]*(1.0-yweight);
                    P_Y[j][yinx+1]  += EwGQ[jj]*(yweight);
                
                }
            }
        }

        // routine to trunc the negligible weights in P_Y, trunc=1e-06
        reduction_matrix();

        // normalize to 1.0.
        sum_PY = 0;
        inv_distri(I_Y, P_Y);
        for(j = 0; j < length_y; j++){sum_PY += I_Y[j]*Y[j];}
        for(j = 0; j < length_y; j++){Y[j]   = Y[j]/sum_PY;}

        printf("\nPermanent process for earnings, states : \n");
        for(y = 0; y<length_y; y++){printf("%f\t",Y[y]);}printf("\n");
        printf("Matrix : \n");
        for(y = 0; y<length_y; y++){for(yy = 0; yy<length_y; yy++){printf("%f\t",P_Y[y][yy]);}printf("\n");}
        printf("\n-- Reduced Matrix : \n");
        for(y = 0; y<length_y; y++){for(yy = 0; yy<max_P_Y[y]; yy++){printf("%d\t",P_Y_reduc[y][yy]);}printf("\n");}
        
    #endif
    
    #if OPT_yprocess == 2

        tauchenfun(rho_y, 1.0, 0.0, std_y, Y, P_Y);    inv_distri(I_Y, P_Y);

        sum_PY = 0.0;
        for(int y = 0; y < length_y; y++){Y[y] = exp(Y[y]);sum_PY += I_Y[y]*Y[y];}
        for(int y = 0; y < length_y; y++){Y[y] = Y[y]/sum_PY;}
        
        printf("\nPersistent process for earnings \n");
        for(y = 0; y<length_y; y++){printf("%f\t",Y[y]);}printf("\n");
        printf("Matrix : \n");
        for(y = 0; y<length_y; y++){for(int y2 = 0; y2<length_y; y2++){printf("%f\t",P_Y[y][y2]);}printf("\n");}printf("\n");

    #endif



    // Transitory process for earnings.
    // This part is taken from Hubmer et al. (2020) -- it combines an unemployment state with a transitory part.
    double E_tmp[length_e-1] = {0}, P_E_tmp[length_e-1] = {0};
    if(length_e == 1){
        E[0]   = 1.0; P_E[0] = 1.0;
    }else{
        // gaussian - quadrature (weights and nodes)
        gauher(E_tmp,P_E_tmp,length_e-1);
        sum_PY = 0;
        for(j = 0; j < length_e-1; j++){sum_PY     += P_E_tmp[j];}
        for(j = 0; j < length_e-1; j++){P_E_tmp[j]  = P_E_tmp[j]/sum_PY;}
        for(j = 0; j < length_e-1; j++){E_tmp[j]    = exp(-pow(std_e,2.0)/2.0 + std_e*E_tmp[j]*sqrt(2.0));}
        
        // unemployment state.
        E[0]    = u_ben;
        P_E[0]  = urate;
        
        // other states.
        for(j = 1; j < length_e; j++){
            E[j]   = E_tmp[j-1];
            P_E[j] = (1.0-P_E[0])*P_E_tmp[j-1];
        }

        // Rescale & Normalisation
        sum_PY = 0;
        for(j = 0; j < length_e; j++){sum_PY += P_E[j];}
        for(j = 0; j < length_e; j++){P_E[j]  = P_E[j]/sum_PY;}
        sum_PY = 0;
        for(j = 0; j < length_e; j++){sum_PY += P_E[j]*E[j];}
        for(j = 0; j < length_e; j++){E[j]    = E[j]/sum_PY;}
    }
    printf("\nTransitory process for earnings \n");
    for(e = 0; e<length_e; e++){printf("%f\t",E[e]);}printf("\n");
    printf("Matrix : \n");
    for(e = 0; e<length_e; e++){printf("%f\t",P_E[e]);}printf("\n");



    /** ------------------ INVESTMENT PRODUCTIVITY ------------------ **/
    // Calibrated to match the distribution of fixed effects over returns in Fagereng (2020, ECMA).
    
    // persistent component.
//    double P_THETA_tmp[length_t][length_t]  = {{1.0-proba_investor,proba_investor,0.0},
//                         {proba_not_investor,1.0-proba_not_investor-prob_move_up,prob_move_up},
//                         {proba_not_investor,prob_move_down,1.0-proba_not_investor-prob_move_down}};
//
//    for(t = 0; t<length_t; t++){for(int t2 = 0; t2<length_t; t2++){P_THETA[t][t2] = P_THETA_tmp[t][t2];}}
                         
    printf("\nPersistent process of investment shock \n");
    for(t = 0; t<length_t; t++){printf("%f\t",THETA[t]);}printf("\n");
    printf("Matrix : \n");
    for(t = 0; t<length_t; t++){for(int t2 = 0; t2<length_t; t2++){printf("%f\t",P_THETA[t][t2]);}printf("\n");}printf("\n");

    // idiosyncratic shock.
    printf("\nIdiosyncratic process of investment shock \n");
    for(z = 0; z<length_z; z++){printf("Z[%d] = %f \n",z,Z[z]);}printf("\n");
    printf("Matrix : \n");
    for(z = 0; z<length_z; z++){printf("%f\t",TZ[z]);}printf("\n");




    /** ------------------ GRID INITIALIZATION ------------------ **/
    for(i=0;i<length_i;i++){PSI[i]       =  gridspace(i,PSI_min,PSI_max,length_i);}     // cash on hand
    for(i=0;i<length_i;i++){A[i]         =  gridspace(i,A_min,A_max,length_i);}         // liquid savings
    for(nu=0;nu<length_nu;nu++){NUgrid[nu]   =  gridspace(nu,0,1.0,length_nu);}//printf("%f \n",NUgrid[nu]);}             // share in equity
    



    /** ------------------ PARAMETERS OF THE MODEL + TEMP FILE ------------------ **/
    // specify the parameters.
    printf("\n");
    for(int p=0; p<nb_para; p++){printf("PARAMETERS:: %s = %1.4f \n", set_params_name[p], set_params[p]);}
    printf("\n\n");

    // specify the moments.
    if(set_obsmoments[0] == 0){printf("error? no observed moments specified for m=0\n\n");getchar();}
    if(set_params[0] == 0){printf("error? no parameters specified p=0\n\n");getchar();}

    // temp file.
    time_t rawtime,timeofstart,timeofend;struct tm * timeinfo;time ( &rawtime );time (&timeofstart);timeinfo = localtime ( &rawtime );
    char buffer [80]; FILE *tempfileoutfile; timeinfo = localtime ( &rawtime );
    strftime (buffer,80,"OUTPUT/temp/tempfile_%Y%m%d@%H%M%S.out",timeinfo);strncpy(FILE_log, buffer, sizeof(FILE_log) - 1);
    FILE_log[sizeof(FILE_log) - 1] = 0;
    tempfileoutfile=fopen(FILE_log, "w");setbuf ( tempfileoutfile , NULL );
        fprintf(tempfileoutfile,"====== Program Starting at %s  \n \n",asctime(timeinfo));
        fprintf(tempfileoutfile,"------ Model parameters ------\n");
        for(int p=0; p<nb_para; p++){fprintf(tempfileoutfile,"%s = %1.4f\n",set_params_name[p],set_params[p]);}
        fprintf(tempfileoutfile,"\n\n");
    fclose(tempfileoutfile);

}





