/** -------------------------------------------------
    Alexandre GAILLARD - PHILIPP WANGNER  (GIF PAPER, 2020)
    
    -------
    POLICY FILE ::
    This file compute the policy functions of the model.
    
    - K_NU_sol      :: compute the optimal share of capital.
    - POLICY_EGM    :: compute the policy functions of the model.
        - cons  :: consumption level
        - save  :: liquid saving
        - NU    :: share invested in risky assets.
        - K     :: investment in equity  (nu*A)
    -------
**/



/// FUNCTION OPTIMAL SHARE OF RISKY ASSET.
struct struct_NU{double *G_VF_struct; int igrid, ygrid, tgrid, NU_stop;};

double FOC_nu(const double nu, void *params){

    struct struct_NU *NUPARAMS   = (struct struct_NU *) params;
    double *G_VF                = NUPARAMS->G_VF_struct;
    int    i                    = NUPARAMS->igrid;
    int    y                    = NUPARAMS->ygrid;
    int    t                    = NUPARAMS->tgrid;

    double G_V_next, weight, psinext, EG_V=0.0;
    int ipsi, z;
    
    // constrain individuals to invest if they generate negative wealth!
    NUPARAMS -> NU_stop = 0.0;
    
    // computing the next period certainty equivalent operator.
    for(z=0;z<length_z;z++){                // investment shock on private equity

        // compute next period wealth.
//        psinext = nu*(A[i]+W*Y[y])*((1.0+riskfree+premium)*Z[z]) + (1.0-nu)*(A[i]+W*Y[y])*(1.0+riskfree);
        psinext = W*Y[y]*ZY[z] + nu*(A[i])*((1.0+riskfree+premium)*ZR[z]) + (1.0-nu)*(A[i])*(1.0+riskfree);
        if(psinext < 0.0){printf("ERROR");}
//        if(psinext > PSI_max){
//            printf("i = %d / psinext = %f,  %f",i, psinext, PSI[i]);getchar();
//            NUPARAMS -> NU_stop=1.0;
//            break;
//        }

        ipsi    = max(0,min((int)(floor(invgridspace(psinext,PSI_min,PSI_max,length_i))),length_i-2));
        weight  = (psinext - PSI[ipsi])/(PSI[ipsi+1] - PSI[ipsi]);
        if(weight > 1){weight = 1.0 + (weight-1.0)/1.2;}    // this is a small correction for potential approx in the extrapolation at the top.
        if(ipsi < 0 | weight < 0.0){printf("mistake: %d %f psinext = %f", ipsi, weight, psinext);getchar();}

        // correct for bounds.
        G_V_next    = inter1d(weight,G_VF[x(ipsi,y,t)],G_VF[x(ipsi+1,y,t)]);   // next period VF.
        EG_V        += TZ[z]*G_V_next;
        
    } // end z
    
    EG_V = (invGCARA(EG_V,A[i],THETA[t]));
    
    return(EG_V);
}







/** Function Policy EGM **/
void POLICY_EGM(double *VF, double *dVF, double *save, double *cons, double *NU, double *X, double *N, int transi_dyn){


    // initialize index
    int    iter,yt,y,t,i,tnext,ynext,threshold_ii,ii,ipsi,nu,nn;
    double CRIT_VF;
    double EVF_next,dEVF_next,psi_next,weight,V_max,tempvf,wS,wV,max_GVF,NUmax,val_try,currentcoh;
    
    
    // File of interest.
    FILE *pfile;


    // initialize pointer
    double *PSI_endo, *EVF_endo, *VF_tilde, *VF_new, *G_VF;
    PSI_endo     =   (double *) calloc((length_x), sizeof(double));       // endogenous grid coh.
    EVF_endo     =   (double *) calloc((length_x), sizeof(double));       // endogenous grid values function next.
    VF_tilde     =   (double *) calloc((length_x), sizeof(double));       // value second sub-period.
    VF_new       =   (double *) calloc((length_x), sizeof(double));       // value new current period.
    G_VF         =   (double *) calloc((length_x), sizeof(double));       // the certainty equivalent operator.

        
    // loop for fixed point (contraction mapping theorem)
    CRIT_VF     =   1.0;
    iter        =   0;

    int stop_transi = 0;
    while ((CRIT_VF > EPS_VALUE) && (iter < 1000) && (stop_transi == 0)){

        if(iter == 999){printf("\n Error VF iteration: maximum iteration reached with critere = %f\n", CRIT_VF);}


        /** ------------------ STEP 1 ------------------
            -   compute the optimal consumption-saving problem given "k" and "psi'"
            -------------------------------------------- **/
            
        // parallelizing
        #if OPT_OMP == 1
            #pragma omp parallel for private (yt,y,t,i,psi_next,ipsi,weight,EVF_next,dEVF_next)
        #endif
        
        for(yt = 0; yt < length_yt; yt++){
        
            y = (int)(floor((yt/(length_t))));       // find index y
            t = (int)(yt-length_t*y);                // find index t
            
            for(i = 0; i < length_i; i++){      // this is the end-of-period decision of saving' (how much to save)
            
                // compute expected value next period (to change with bequest later.)
                EVF_next  = 0;
                dEVF_next = 0;
                for(ynext = 0; ynext < length_y; ynext++){
                    for(tnext = 0; tnext < length_t; tnext++){
                        EVF_next += beta*((1.0-P_die)*P_Y[y][ynext]*P_THETA[t][tnext]*(VF[x(i,ynext,tnext)]));
                        dEVF_next += beta*((1.0-P_die)*P_Y[y][ynext]*P_THETA[t][tnext]*(dVF[x(i,ynext,tnext)]));
                    }  // end tnext
                } // end ynext
                
                // compute the current cash-on-hand.
                PSI_endo[x(i,y,t)] = A[i] + invMU(dEVF_next);  // unconstrained problem :: psi = a' + c
                EVF_endo[x(i,y,t)] = EVF_next;
                if(PSI_endo[x(i,y,t)] != PSI_endo[x(i,y,t)]){printf("bug COH : %d %d %d | %f %f %f %f \n ", i,y,t,PSI_endo[x(i,y,t)], invMU(EVF_next),EVF_endo[x(i,y,t)],A[i]); getchar();}
            
            }   // end i
        }   // end yt
        

        // print in log-files.
        #if PRINT_SS == 1
//            printf("\nWrite endogenous grid into log-file \n");
            pfile=fopen(FILE_endopolicy, "w"); setbuf(pfile, NULL );
            for(y=0;y<length_y;y++){
                for(t=0;t<length_t;t++){
                    for(i=0;i<length_i;i++){
                        fprintf(pfile,"%5d\t%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\n",i,t,y,A[i],PSI_endo[x(i,y,t)],EVF_endo[x(i,y,t)]);
                    }
                }
            }
            fclose(pfile);
        #endif



        /** ------------------ STEP 2 ------------------
            -   Rescale the solution.
            -------------------------------------------- **/
        
        // parallelizing
        #if OPT_OMP == 1
            #pragma omp parallel for private (yt,y,t,threshold_ii,i,tempvf,wV,wS,ii,weight,currentcoh)
        #endif
        
        for(yt = 0; yt < length_yt; yt++){
    
            y = (int)(floor((yt/(length_t))));       // find index y
            t = (int)(yt-length_t*y);                // find index t
        
            threshold_ii = 0;
            for(i = 0; i < length_i; i++){  // current grid of A.

                currentcoh = PSI[i];

                // case 1: eat the borrowing limit.
                if(currentcoh < PSI_endo[x(0,y,t)]){save[x(i,y,t)] = A[0];  tempvf = EVF_endo[x(0,y,t)];}
                
                // case 2: extrapolation.
                if(currentcoh >= PSI_endo[x(length_i-1,y,t)]){
                    wV  = (EVF_endo[x((length_i-1),y,t)]-EVF_endo[x((length_i-2),y,t)])/(PSI_endo[x(length_i-1,y,t)]-PSI_endo[x(length_i-2,y,t)]);
                    wS  = (A[(length_i-1)]-A[(length_i-2)])/(PSI_endo[x(length_i-1,y,t)]-PSI_endo[x(length_i-2,y,t)]);
                    
                    save[x(i,y,t)]  = (currentcoh-PSI_endo[x(length_i-1,y,t)])*wS + A[(length_i-1)];
                    tempvf          = (currentcoh-PSI_endo[x(length_i-1,y,t)])*wV + EVF_endo[x((length_i-1),y,t)];
                }
                
                // case 3: interior solution.
                if(currentcoh < PSI_endo[x(length_i-1,y,t)] && currentcoh >= PSI_endo[x(0,y,t)]){
                
                    ii      = max(threshold_ii,0);
                    
                    while((currentcoh>PSI_endo[x(ii,y,t)]) && (ii < length_i-1)){ii++;}
                    threshold_ii = max(ii-2,0);
                    
                    if(ii == length_i-1 && currentcoh>PSI_endo[x(ii,y,t)]){printf("mistake 4: %f %f %f %d %d \n", currentcoh, PSI_endo[x(length_i-2,y,t)], PSI_endo[x(length_i-1,y,t)], ii, i);getchar();}
                
                    weight           = (currentcoh - PSI_endo[x((ii-1),y,t)])/(PSI_endo[x(ii,y,t)] - PSI_endo[x((ii-1),y,t)]);
                    save[x(i,y,t)]   = inter1d(weight,A[(ii-1)],A[ii]);
                    tempvf           = inter1d(weight,EVF_endo[x((ii-1),y,t)],EVF_endo[x(ii,y,t)]);
                }
                
                cons[x(i,y,t)]       = currentcoh - save[x(i,y,t)] - A_min;
                VF_tilde[x(i,y,t)]   = U(cons[x(i,y,t)]) + tempvf;
                G_VF[x(i,y,t)] = GCARA((VF_tilde[x(i,y,t)]),currentcoh,THETA[t]);  // extra value in case we want to make risk aversion type/wealth dep
                
            }
        }
    
        
        // print in log-files.
        #if PRINT_SS == 1
//            printf("\nWrite policy function (second sub-period) into log-file \n");
            pfile   =   fopen(FILE_policy2, "w"); setbuf(pfile, NULL );
            for(yt = 0; yt < length_yt; yt++){
                y = (int)(floor((yt/(length_t))));       // find index y
                t = (int)(yt-length_t*y);                // find index t
                for(i=0;i<length_i;i++){
                    fprintf(pfile,"%5d\t%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",i,t,y,currentcoh,(invGCARA(G_VF[x(i,y,t)],currentcoh,THETA[t])), G_VF[x(i,y,t)],VF_tilde[x(i,y,t)],cons[x(i,y,t)],save[x(i,y,t)]);
                }
            }
            fclose(pfile);
        #endif
        
        
            
            
        /** ------------------ STEP 3 ------------------
            -   compute the optimal K and share of risky assets.
            -   given the second sub-period solution VF_tilde, cons, save, solve for the optimal k, risky share.
            ---------------------------------------------------------------------------------------------------------- **/
            
        struct struct_NU params_NU;
        CRIT_VF = 0.0;
        
        // parallelizing
        #if OPT_OMP == 1
            #pragma omp parallel for reduction(max:CRIT_VF) private (params_NU,nu,nn,yt,y,t,i,max_GVF,NUmax,val_try)
        #endif
        
        for(yt = 0; yt < length_yt; yt++){
        
            y = (int)(floor((yt/(length_t))));       // find index y
            t = (int)(yt-length_t*y);                // find index t
     
            for(i = 0; i < length_i; i++){  // next period grid of A.
            
                // --- OPTIMAL INVESTMENT SHARE
                // set the params.
                params_NU.G_VF_struct      = G_VF;
                params_NU.ygrid            = y;
                params_NU.igrid            = i;
                params_NU.tgrid            = t;
                
                
                // compute the optimal share invested.
                max_GVF = -100000;
                for(nn = 0; nn<length_nu; nn++){
                    val_try = FOC_nu(NUgrid[nn],&params_NU);

                    if(val_try > max_GVF){
                        max_GVF = val_try;
                        NUmax   = NUgrid[nn];
                    }
                }

                VF_new[x(i,y,t)] = max_GVF;
                NU[x(i,y,t)]     = NUmax;

//                    NU[x(i,y,t)]     = 0.20;
//                    VF_new[x(i,y,t)] = FOC_nu(NU[x(i,y,t)],&params_NU);

                // compute the criterion
                if(i < length_i){CRIT_VF = max(CRIT_VF,fabs(VF[x(i,y,t)] - VF_new[x(i,y,t)]));}
                    
                // replace the value (contraction mappin apply)
                VF[x(i,y,t)] = VF_new[x(i,y,t)];
                
                // compute the derivative.
                derivative(i,y,t,VF,dVF);

            } // end i
        } // end yt
            


        // print in log-files.
        #if PRINT_SS == 1
            pfile   =   fopen(FILE_policy, "w"); setbuf(pfile, NULL );
            for(yt = 0; yt < length_yt; yt++){
                y = (int)(floor((yt/(length_t))));       // find index y
                t = (int)(yt-length_t*y);                // find index t
                for(i=0;i<length_i;i++){
                    fprintf(pfile,"%5d\t%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\n",i,t,y,A[i],VF[x(i,y,t)],NU[x(i,y,t)]);
                }
            }
            fclose(pfile);
        #endif
    
    
        printf("ITERATION = %d,  CRIT_VF = %20.15f\n", iter, CRIT_VF); //getchar();

        if(transi_dyn == 1){stop_transi = 1;}
        iter++;
        
    }//end of while loop






    // print in log-files.
    #if PRINT_SS == 2
        pfile=fopen(FILE_endopolicy, "w"); setbuf(pfile, NULL );
        for(y=0;y<length_y;y++){
            for(t=0;t<length_t;t++){
                for(i=0;i<length_i;i++){
                    fprintf(pfile,"%5d\t%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\n",i,t,y,A[i],PSI_endo[x(i,y,t)],EVF_endo[x(i,y,t)]);
                }
            }
        }
        fclose(pfile);
        
        
        
        pfile   =   fopen(FILE_policy2, "w"); setbuf(pfile, NULL );
        for(yt = 0; yt < length_yt; yt++){
            y = (int)(floor((yt/(length_t))));       // find index y
            t = (int)(yt-length_t*y);                // find index t
            for(i=0;i<length_i;i++){
                fprintf(pfile,"%5d\t%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",i,t,y,A[i],VF_tilde[x(i,y,t)],cons[x(i,y,t)],save[x(i,y,t)]);
            }
        }
        fclose(pfile);
        
        
        
        printf("\n Write policy function into log-file \n");
        pfile   =   fopen(FILE_policy, "w"); setbuf(pfile, NULL );
        for(yt = 0; yt < length_yt; yt++){
            y = (int)(floor((yt/(length_t))));       // find index y
            t = (int)(yt-length_t*y);                // find index t
            for(i=0;i<length_i;i++){
                fprintf(pfile,"%5d\t%5d\t%5d -- \t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",i,t,y,A[i],VF[x(i,y,t)],NU[x(i,y,t)],(save[x(i,y,t)]));
            }
        }
        fclose(pfile);
    #endif

    getchar();



    // free up the pointers.
    free(PSI_endo);
    free(VF_new);
    free(VF_tilde);
    free(EVF_endo);
    free(G_VF);

}


