/* ********************************************** */
//      CONTROL RANDOM SEARCH -- VERSION 2        //
/*                    JUNE 2018                   */
//               Alexandre Gaillard               //
/* ********************************************** */

// THIS CAN BE APPLIED ON CLUSTER SERVER USING DIFFERENT NODE //
// Nodes successively replace the worst point in the set of known parameters //
// converge to the global minimum //



int rand_a_b(int a, int b){
    return rand()%(b-a) + a;
}



bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

using namespace std;





double myCRS(double *best_point, double fun(double *)) {
    
    
    /******************************/
    /**   GENERATE A SOBOL SEQ   **/
    /******************************/
    
    int nbNODE, s, p;
    cout << "NODE number (integer): ";
    cin >> nbNODE;
    
    printf("\n \n Number of the node: %d \n \n ", nbNODE);
    
    FILE *soboloutfile;

    /******** INITIALIZE SOBOL SEQUENCE ********/
    char buffer [80];
    int numbersobol = 0;
    sprintf(buffer,"INPUT/SMM_para/SOBOL_%d_NODE_%d.out",numbersobol,nbNODE);
    strncpy(SOBOL, buffer, sizeof(SOBOL) - 1);
    SOBOL[sizeof(SOBOL) - 1] = 0;
    

    if(is_file_exist(SOBOL)){
    
        printf("already exit \n \n ");
        // CAUTION, dimension must be placed in the right order of the preceding loop constructing the SOBOL sequence //
        // here it is: dim1 = nbsovol, dim2 = nb_para
        readINPUT2D(SOBOL_seq, SOBOL);
    
    }else{
    #if SOBOLCREATE == 1
        printf("do not already exit \n \n ");
        
        //SOBOL GENERATION
        float *sobolinitial;
        sobolinitial = i4_sobol_generate(nbsobol, nb_para, (50*(nbNODE+1)));

        //SAVE Sobol
        soboloutfile=fopen(SOBOL, "w");
        setbuf ( soboloutfile , NULL );

        for(s = 0; s < nbsobol; s++){
            for(p = 0; p < nb_para; p++){
                SOBOL_seq[s][p] = sobolinitial[s*nb_para + p];
                fprintf(soboloutfile,"%20.15f\n", SOBOL_seq[s][p]);
            }
        }
        
        fclose(soboloutfile);
    #endif
    }
    
    
    
    
    
    
    
    /*******************************/
    /**   CONTROL RANDOM SEARCH   **/
    /*******************************/
    
    // CRS2 - as in Kaelo and Ali (2006) "Some variants of the controlled random search algorithm for global optimization"


    // initialization of variables
    bool STOP = false;
    bool belong_set;
    double new_point[nb_para];
    double new_coordinate;
    double new_evaluation;
    
    double POINTS_temp[nb_para];
    double EVALUATION[nbsobol];             // evaluation of the function at the N points
    
    int index_low, index_high, index;

    
    // for sample generation:
    int SAMPLE[nb_para+1];
    std::uniform_int_distribution<> unif_dist(0, (nbsobol-1));
    std::mt19937 gen(123*nbNODE); //
    
    
    

    
    
    // Step 0: Initialization
    // ----------------------
 
    /** 0. REWRITE THE SOBOL_SEQ IN TERMS OF PARAMETERS **/

    // b) Compute the corresponding points from the sequence
    for(s = 0; s < nbsobol; s++) {
        for(p = 0; p < nb_para; p++) {
            PARA_seq[s][p] = set_params_min[p] + SOBOL_seq[s][p]*(set_params_max[p] - set_params_min[p]);
        }
    }


    for(s = 0; s < nbsobol; s++){
        printf("SOBOL: %d \t", s);
    }
    printf("\n");
    for(p = 0; p < nb_para; p++) {
        for(s = 0; s < nbsobol; s++) {
            printf("%f \t", PARA_seq[s][p]);
        }
        printf("\n");
    }
    


    
    
    /** 0.1 Evaluate the function at the N generated points **/
    for(s = 0; s < nbsobol; s++) {
    
        printf("\n \n SOBOL NUMBER:: %d \n \n", s);
        
        for(p = 0; p < nb_para; p++) {
            POINTS_temp[p] = PARA_seq[s][p];
        }
        
        EVALUATION[s] = fun(POINTS_temp);    // compute the equilibrium for each NODE
    }
    
    

    // write the new sobol //
    char evalwrite[200],paramswrite[200], replacewrite[200];
    sprintf(evalwrite,"INPUT/node/evaluation_initial_node_%d.out",nbNODE);
    strncpy(evalwrite, evalwrite, sizeof(evalwrite) - 1);
    evalwrite[sizeof(evalwrite) - 1] = 0;
    
    sprintf(paramswrite,"INPUT/node/params_initial_node_%d.out",nbNODE);
    strncpy(paramswrite, paramswrite, sizeof(evalwrite) - 1);
    paramswrite[sizeof(paramswrite) - 1] = 0;

    FILE *EVALfile, *PARAMSfile;
    EVALfile=fopen(evalwrite, "w");setbuf ( EVALfile , NULL );
    for(s = 0; s < nbsobol; s++) {
        fprintf(EVALfile,"%20.15f\n", EVALUATION[s]);
    }
    fclose(EVALfile);
    
    PARAMSfile=fopen(paramswrite, "w");setbuf ( PARAMSfile , NULL );
    for(s = 0; s < nbsobol; s++) {
        for(p = 0; p < nb_para; p++) {
           fprintf(PARAMSfile,"%20.15f\n", PARA_seq[s][p]);
        }
    }
    fclose(PARAMSfile);
    
    

 
 
    // FIRST STEP OF REPLACEMENT OF THE BEST POINTS //
    // All nodes that end here should compare their EVALUATION set to the existing ONE //
    // A node has to compare its sample S with the best sample hat(S) //
    // File evaluation_set.out  contains the best evaluation of all the nodes //
    // File parameter_set.out contains the best points of all the nodes //
    
    /** 1. READ THE EXISTING SAMPLE hat(S) IF EXIST **/
    if(is_file_exist("INPUT/CRS/evaluation_set.out")){
      
        double EVALUATION_compare[nbsobol],PARA_compare[nbsobol][nb_para], REPLACEMENT;
        int max_dist;
      
      
        readINPUT1D(EVALUATION_compare,"INPUT/CRS/evaluation_set.out");
        readINPUT2D(PARA_compare,"INPUT/CRS/parameter_set.out");
      
      
        // CHECK THE FILES //
        for(s = 0; s < nbsobol; s++) { // eval by eval
    
            max_dist = 0;
            REPLACEMENT = 0;
            for(int n = 0; n < nbsobol; n++) { // we want to replace the less efficient one !
                if(EVALUATION[s] < EVALUATION_compare[n] && EVALUATION_compare[n] >= EVALUATION_compare[max_dist]){  // check if eval is better
                    max_dist = n;
                    REPLACEMENT = 1;
                }
            } // end check if better point in the set.
            
            if(REPLACEMENT == 1){
                EVALUATION_compare[max_dist] = EVALUATION[s];   // replace the point "max_dist"
                for(p = 0; p < nb_para; p++) {
                    PARA_compare[max_dist][p] = PARA_seq[s][p];    // replace the points in set "n"
                }
            }
            
        } // end check eval
        
        
        // REWRITE THE NEW FILES //
        char file_name[200];
        FILE *write_file;
        
        // evaluation_set.out
        sprintf(file_name,"INPUT/CRS/evaluation_set.out");
        strncpy(file_name, file_name, sizeof(file_name) - 1); file_name[sizeof(file_name) - 1] = 0;

        write_file=fopen(file_name, "w"); setbuf ( write_file , NULL );

        for(s = 0; s < nbsobol; s++) {
            fprintf(write_file,"%20.15f\n",EVALUATION_compare[s]);
        }
        
        fclose(write_file);
        
        
        // parameter_set.out
        sprintf(file_name,"INPUT/CRS/parameter_set.out");
        strncpy(file_name, file_name, sizeof(file_name) - 1); file_name[sizeof(file_name) - 1] = 0;

        write_file=fopen(file_name, "w"); setbuf ( write_file , NULL );

        for(s = 0; s < nbsobol; s++) {
            for(p = 0; p < nb_para; p++) {
                fprintf(write_file,"%20.15f\n",PARA_compare[s][p]);
            }
        }
    
        fclose(write_file);
        
        
        // UPDATE THE SEQUENCE EVALUATION AND PARAMETERS //
        for(s = 0; s < nbsobol; s++) {
            for(p = 0; p < nb_para; p++) {
                PARA_seq[s][p] = PARA_compare[s][p];
            }
            
            EVALUATION[s] = EVALUATION_compare[s];    // compute the equilibrium for each NODE
        }
        
        
    }else{
 
        // WRITE THE NEW FILES //
        char file_name[200];
        FILE *write_file;
        
        // evaluation_set.out
        sprintf(file_name,"INPUT/CRS/evaluation_set.out");
        strncpy(file_name, file_name, sizeof(file_name) - 1); file_name[sizeof(file_name) - 1] = 0;

        write_file=fopen(file_name, "w"); setbuf ( write_file , NULL );

        for(s = 0; s < nbsobol; s++) {
            fprintf(write_file,"%20.15f\n",EVALUATION[s]);
        }
        
        fclose(write_file);
        
        
        // parameter_set.out
        sprintf(file_name,"INPUT/CRS/parameter_set.out");
        strncpy(file_name, file_name, sizeof(file_name) - 1); file_name[sizeof(file_name) - 1] = 0;

        write_file=fopen(file_name, "w"); setbuf ( write_file , NULL );

        for(s = 0; s < nbsobol; s++) {
            for(p = 0; p < nb_para; p++) {
                fprintf(write_file,"%20.15f\n",PARA_seq[s][p]);
            }
        }
    
        fclose(write_file);
    
    } // file exist
    
    
    
    
    
    
    // BEGIN LOOP  (this will be done for each NODE) //
    int iteration_count = 0;
    while(STOP == false) {

        /** Step 1: Find the best and worst points (within the sequence of the NODE) **/
        // --------------------------------------
        // Stop if f_h - f_l < epsilon


        index_low = 0; index_high = 0; // to initialize the first time
        
        for(s = 0; s < nbsobol; s++) {
            if(EVALUATION[s] <= EVALUATION[index_low]) { index_low = s; } // EVALUATE[index_low] = f_l
            if(EVALUATION[s] >= EVALUATION[index_high]) { index_high = s; } // EVALUATE[index_low] = f_h
        }
        

        /** Step 2: Stopping rule **/
        // ---------------------
        
        if((EVALUATION[index_high] - EVALUATION[index_low]) < criterionCRS) { STOP = true; }

        printf("\n \n CRS DISTANCE: == %f \n \n ",  (EVALUATION[index_high] - EVALUATION[index_low]));


        /** Step 3: Generate a trial point (randomly drawn) **/
        // ------------------------------

        belong_set = false; // to initialize, say that the point does not belong to the set.

        while(belong_set == false) { // as long as the generated trial point does not belong to the set, do NOT evaluate the function

            //---- 3.1 Randomly pick n+1 points ----//
            // ---------------------------
            // remark: CRS2 method, so always include the set of points with the lowest value in this set (as the first point).
            // just pick n+1 index among the set 1:N (i.e. 0:(N-1) on C++)
            // <=> Generate n+1 different index (for CRS) (in reality only N-1 index generated the first one being the lowest index)

            SAMPLE[0] = index_low; // always include it as first point (CRS2)
            


            // Generate a contraction using nb_para+1 (n+1) << nbsobol (N) with replacement (randomly picked) //
            for(int n = 1; n < (nb_para + 1); n++){
                SAMPLE[n] = unif_dist(gen);       // will be different for each NODE due to "gen"
            }



            //---- 3.2 Generate new trial points (that will replace the worst point later) ----//
            for(p = 0; p < nb_para; p++) { // param by param
                new_coordinate = 0;
                for(int i=0; i < (nb_para + 1); i++) { // point in sample by point in sample  (we choose sample = nb_para + 1)
                    if(i != (nb_para)) {
                        new_coordinate += (2*PARA_seq[SAMPLE[i]][p])/(nb_para); //  the first points: 2 * mean (centroid)
                    } else {
                        new_coordinate = new_coordinate - PARA_seq[SAMPLE[i]][p]; // for the last point we substract its value
                    }
                }
                new_point[p] = new_coordinate;
            }

            
            // ---- 3.3 Check if the point is within our Set limits ---- //
            // if does not belong to set, will just try some other point //
            
            belong_set = true;
            for(p = 0; p < nb_para; p++) {
                if(new_point[p] < set_params_min[p]) { belong_set = false; } // does not belong to set in this case
                if(new_point[p] > set_params_max[p]) { belong_set = false; }
            }

   
            // End-up with a new set of parameters that belongs to the initial sets //
            // For now, we did not excluded inefficient points //
            
        } // end while on belong_set
  
  

        /** 4. Evaluate function at this new point **/
        new_evaluation = fun(new_point);
        
        
        
        /** 5. UPDATE THE KNOWLEDGE ABOUT BEST POINTS **/
        double EVALUATION_compare[nbsobol],PARA_compare[nbsobol][nb_para], REPLACEMENT;
        int max_dist;
      
        readINPUT1D(EVALUATION_compare,"INPUT/CRS/evaluation_set.out");
        readINPUT2D(PARA_compare,"INPUT/CRS/parameter_set.out");
        
        
        // CHECK THE FILES //
        max_dist = 0;
        REPLACEMENT = 0;
        for(s = 0; s < nbsobol; s++) { // we want to replace the less efficient one !
            if(new_evaluation < EVALUATION_compare[s] && EVALUATION_compare[s] >= EVALUATION_compare[max_dist]){  // check if eval is better
                max_dist = s;
                REPLACEMENT = 1;
            }
        } // end check if better point in the set.
        
        
        
        /** 5. Replace if the evaluation is below the highest one **/
        if(REPLACEMENT == 1){
        
            EVALUATION_compare[max_dist] = new_evaluation;   // replace the point "max_dist"
            for(p = 0; p < nb_para; p++) {
                PARA_compare[max_dist][p] = new_point[p];    // replace the points in set "n"
            }

            
            // REWRITE THE NEW FILES //
            char file_name[200];
            FILE *write_file;
            
            // evaluation_set.out
            sprintf(file_name,"INPUT/CRS/evaluation_set.out");
            strncpy(file_name, file_name, sizeof(file_name) - 1); file_name[sizeof(file_name) - 1] = 0;

            write_file=fopen(file_name, "w"); setbuf ( write_file , NULL );

            for(s = 0; s < nbsobol; s++) {
                fprintf(write_file,"%20.15f\n",EVALUATION_compare[s]);
            }
            
            fclose(write_file);
            
            
            // parameter_set.out
            sprintf(file_name,"INPUT/CRS/parameter_set.out");
            strncpy(file_name, file_name, sizeof(file_name) - 1); file_name[sizeof(file_name) - 1] = 0;

            write_file=fopen(file_name, "w"); setbuf ( write_file , NULL );

            for(s = 0; s < nbsobol; s++) {
                for(p = 0; p < nb_para; p++) {
                    fprintf(write_file,"%20.15f\n",PARA_compare[s][p]);
                }
            }
        
            fclose(write_file);
        
            
            // UPDATE THE SEQUENCE EVALUATION AND PARAMETERS //
            for(s = 0; s < nbsobol; s++) {
            
                for(p = 0; p < nb_para; p++) {
                    PARA_seq[s][p] = PARA_compare[s][p];
                }
                
                EVALUATION[s] = EVALUATION_compare[s];    // compute the equilibrium for each NODE
            }
        
        } // otherwise no replacement
        
 

        // temp file to save where I am //
        sprintf(replacewrite,"INPUT/node/replacement_node=%d.out",nbNODE);
        strncpy(replacewrite, replacewrite, sizeof(replacewrite) - 1);
        replacewrite[sizeof(replacewrite) - 1] = 0;
        
        sprintf(evalwrite,"INPUT/node/evaluation_node=%d_%d_highest=%d_lowest=%d.out",nbNODE,iteration_count,index_high,index_low);
        strncpy(evalwrite, evalwrite, sizeof(evalwrite) - 1);
        evalwrite[sizeof(evalwrite) - 1] = 0;
        
        sprintf(paramswrite,"INPUT/node/params_node=%d_%d.out",nbNODE,iteration_count);
        strncpy(paramswrite, paramswrite, sizeof(evalwrite) - 1);
        paramswrite[sizeof(paramswrite) - 1] = 0;

        FILE *REPLACEfile;
        REPLACEfile=fopen(replacewrite, "a");setbuf ( REPLACEfile , NULL );
        fprintf(REPLACEfile,"%20.15f\n", REPLACEMENT);
        fclose(REPLACEfile);
        
        EVALfile=fopen(evalwrite, "w");setbuf ( EVALfile , NULL );
        for(s = 0; s < nbsobol; s++) {
            fprintf(EVALfile,"%20.15f\n", EVALUATION[s]);
        }
        fclose(EVALfile);
        
        PARAMSfile=fopen(paramswrite, "w");setbuf ( PARAMSfile , NULL );
        for(s = 0; s < nbsobol; s++) {
            for(p = 0; p < nb_para; p++) {
               fprintf(PARAMSfile,"%20.15f\n", PARA_seq[s][p]);
            }
        }
        fclose(PARAMSfile);
        
        


        /** update iteration count **/
        iteration_count = iteration_count + 1; // remark: consider that each function evaluation, even if > max, is a new iteration

        // and then go try a new point
    }
    
    
    for(p = 0; p < nb_para; p++) {
        best_point[p] = PARA_seq[index_low][p];
    }
    
    return EVALUATION[index_low];
    
}
