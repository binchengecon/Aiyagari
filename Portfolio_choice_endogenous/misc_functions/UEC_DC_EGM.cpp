
///***************************************/
///** UPPER ENVELOPE CONDITION (DC-EGM) **/
///***************************************/


// WARNING //

// INTERPOLATE FOR CONSUMPTION, BECAUSE SAVING BEHAVIOR MIGHT CHANGE //
// VERIFICATION FOR WEIGHT, IN SEGMENT 2 //


void UEC_h(double *COHendo, double *Vendo, int *kink_point){

double slope, slope2, coh_temp, coh_temp_next, weight, V_temp;
double slope_J1, slope_J3, coh_kink;
int y, ilim, i, ii, iii, z, t, kk,k, k_J1, k_J3, kkk, i_m, i_p, i_low, i_up, verifcase, i_low_up, Icase, IIcase;



/* Algorithm 3 - Upper Envelope in Iskhakov */

for(y=0;y<length_y;y++){ // current y
    for(h=0; h<length_h; h++){ // current value of h
    
            clim = 0;
            c = 1;
            
            while(c<(length_c_ev)){ // previous k

            
            if(c >= (length_c_ev)){printf("STOP");getchar();}
            if(c <= 0){printf("STOP");getchar();}
                
            verifcase = 0;
            Ccase = 0;
            
            /** 1. DETECT NON CONCAVE REGION -- BREAK THE MONOTONICITY ASSUMPTION **/
            /* Proceed: check whether or not cohendo is increasing everywhere (multiple solutions checking) */
            
            if(COHendo[inxCHY_ev(c,h,y)] < COHendo[inxCHY_ev((c-1),h,y)] && COHendo[inx(i,y,a,z,t)] != -10000){
                //printf("KINK EXIST: %f %f %d %d", COHendo[inx((i-1),y,a,z,t)], COHendo[inx(i,y,a,z,t)], i, a); getchar();
                Icase = 1;
            }
            
            //}

            /** 2. find the point where COHendo(k) < COHendo(k+1) **/
            if(Icase == 1){
                
                ii = i;
                
     
                // Bound the set J2, where there are several value of COH around the same place.
                while((COHendo[inx(ii,y,a,z,t)] > COHendo[inx((ii+1),y,a,z,t)]) && (ii < length_i)){
                    if(ii >= (length_i-2)){IIcase = 2;
                        printf("stop problem in upper envelope condition -- bound reached: %d, %f %f, %d", i, COHendo[inx(ii,y,a,z,t)], COHendo[inx((ii+1),y,a,z,t)], ii);getchar();break;
                    }else{ii++;}
                }
                
//                if(COHendo[inx(i,y,a,z,t)] < COHendo[inx((i-1),y,a,z,t)] && i >= ilim){
//                    printf("BECOME INCREASING FOR: %d, %f %f, %d %f %f", ii, COHendo[inx((ii),y,a,z,t)], COHendo[inx((ii+1),y,a,z,t)], i, COHendo[inx((i-1),y,a,z,t)], COHendo[inx(i,y,a,z,t)]); getchar();
//                }
//
//                if(COHendo[inx(i,y,a,z,t)] < COHendo[inx((i-1),y,a,z,t)] && i >= ilim){
//                    printf("CHECK COH: %f %f %f %f %f %f %f %f %f %f %f \n", COHendo[inx((i-5),y,a,z,t)], COHendo[inx((i-4),y,a,z,t)],COHendo[inx((i-3),y,a,z,t)],COHendo[inx((i-2),y,a,z,t)],COHendo[inx((i-1),y,a,z,t)],COHendo[inx((i),y,a,z,t)],COHendo[inx((i+1),y,a,z,t)],COHendo[inx((i+2),y,a,z,t)],COHendo[inx((i+3),y,a,z,t)],COHendo[inx((i+4),y,a,z,t)],COHendo[inx((i+5),y,a,z,t)]);
//                    printf("CHECK V: %f %f %f %f %f %f %f %f %f %f %f \n", Vendo[inx((i-5),y,a,z,t)], Vendo[inx((i-4),y,a,z,t)],Vendo[inx((i-3),y,a,z,t)],Vendo[inx((i-2),y,a,z,t)],Vendo[inx((i-1),y,a,z,t)],Vendo[inx((i),y,a,z,t)],Vendo[inx((i+1),y,a,z,t)],Vendo[inx((i+2),y,a,z,t)],Vendo[inx((i+3),y,a,z,t)],Vendo[inx((i+4),y,a,z,t)],Vendo[inx((i+5),y,a,z,t)]); getchar();
//                }
//
                
                /** 3. FIND THE REGION "M" where COH cross each other */
                /* FIRST, find left bound */
                i_p = ii;
                i_m = (i-1);
                if(i_m < 0){printf("ERROR 1 %d", i_m);getchar();}
 
                coh_temp_next = COHendo[inx(i_p,y,a,z,t)];
                coh_temp = COHendo[inx(i_m,y,a,z,t)];

                while(coh_temp_next < coh_temp){
                    i_m = i_m - 1;
                    if(i_m == 0){break;}
                    if(i_m < 0){printf("ERROR 2 %d %d %d %f %f, %f %f %f %f %f %f", i_m, i, (i-1), coh_temp_next, coh_temp, COHendo[inx(0,y,a,z,t)], COHendo[inx(1,y,a,z,t)], COHendo[inx(2,y,a,z,t)], COHendo[inx(3,y,a,z,t)], COHendo[inx(4,y,a,z,t)], COHendo[inx(ii,y,a,z,t)]);getchar();}
                    coh_temp = COHendo[inx(i_m,y,a,z,t)];
                }
                
                i_low = i_m;
                
                //if(a == 43){printf("HERE 555 %d %d %d %d", ii, i, i_low, i_up);getchar();}
                
                /* FIRST, find right bound */
                i_p = ii;
                i_m = (i-1);
                if(i_m < 0){printf("ERROR 3 %d", i_m);getchar();}
                
                coh_temp_next = COHendo[inx(i_p,y,a,z,t)];
                coh_temp = COHendo[inx(i_m,y,a,z,t)];
                
                while(coh_temp_next < coh_temp){
                    i_p = i_p + 1;
                    if(i_p > (length_i-1)){printf("ERROR 4 %d", i_p);getchar();}
                    coh_temp_next = COHendo[inx(i_p,y,a,z,t)];
                }
                
                i_up = i_p;
                
                
                // CHECK IF THERE IS NO TWO KINKS OR MORE //
                i_p = i_up;
                i_m = (i_up-1);
                if(i_m < 0){printf("ERROR 5 %d ", i_m);getchar();}
                
                coh_temp_next = COHendo[inx(i_p,y,a,z,t)];
                coh_temp = COHendo[inx(i_m,y,a,z,t)];

                while(coh_temp < coh_temp_next){
                    i_p = i_m;
                    i_m = i_m - 1;
                    if(i_m < 0){printf("ERROR  6 %d", i_m);getchar();}
                    coh_temp_next = COHendo[inx(i_p,y,a,z,t)];
                    coh_temp = COHendo[inx(i_m,y,a,z,t)];
                }
                
                i_low_up = i_p;
                
                
                //if(a == 43){printf("HERE 22 %d %d %d %d", ii, i, i_low, i_up);getchar();}
//                printf("BOUND of M = %d and %d", i_low, i_up);getchar();
 
                
                if(i_low_up != ii){ // case where there are several regions. -- IN that case, INTERPOLATE EVERYTHING ON ONLY TWO SEGMENTS, AND CHECK AFTERWARD IF THIS IS OK.
                    
                    // first, check if there is another lowest value in that segment //
//                    if(COHendo[inx(ii,y,a,z,t)] > COHendo[inx(i_low_up,y,a,z,t)]){i_low_low = (i_low_up);}
//                    if(COHendo[inx(ii,y,a,z,t)] < COHendo[inx(i_low_up,y,a,z,t)]){i_low_low = (ii);}
//
//                    for(k = (ii+1); k<(i_low_up); k++){
//                        if(COHendo[inx(k,y,a,z,t)] < COHendo[inx(i_low_low,y,a,z,t)]){i_low_low = k;}
//                    }
//
//                    // RECHECK the lowest value at the left;
//                    i_p = i_low_low;
//                    i_m = (i-1);
//
//                    coh_temp_next = COHendo[inx(i_p,y,a,z,t)];
//                    coh_temp = COHendo[inx(i_m,y,a,z,t)];
//
//                    while(coh_temp_next < coh_temp){
//                        i_m = i_m - 1;
//                        coh_temp = COHendo[inx(i_m,y,a,z,t)];
//                    }
//
//                    i_low = i_m;
                    
                    
                    
                    // CREATE INTERPOLATED VALUE ON SEGMENT J3 AND J1 ONLY //
                    double V_J1[(i_up - i_low + 1)], V_J3[(i_up - i_low + 1)];
                    
                    
                    for(k = (i_low); k<(i_up+1); k++){
                    
                        // CASE WHERE IN REGION 1 //
                        if(k >= i_low && k <= (i-1)){
                            
                            // VALUE ON SEGMENT 1 IS KNOWN //
                            V_J1[(k-i_low)] = Vendo[inx(k,y,a,z,t)];
                            
                            // VALUE ON SEGMENT 3 IS INTERPOLATED //
                            if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_up),y,a,z,t)]){
                                printf("impossible error 134");getchar();
                            }
                            if(COHendo[inx(k,y,a,z,t)] < COHendo[inx((i_low_up),y,a,z,t)]){
                                slope = (Vendo[inx((i_low_up+1),y,a,z,t)] - Vendo[inx((i_low_up),y,a,z,t)])/(COHendo[inx((i_low_up+1),y,a,z,t)] - COHendo[inx((i_low_up),y,a,z,t)]);
                                V_J3[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i_low_up),y,a,z,t)]) + Vendo[inx((i_low_up),y,a,z,t)];
                            }
                            if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_low_up),y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i_up),y,a,z,t)]){ // INTERPOLATED
                                kk = (i_low_up+1);
                                while(COHendo[inx(k,y,a,z,t)] > COHendo[inx(kk,y,a,z,t)]){
                                    if(kk == (0)){printf("mistake, impossible case 0145");getchar();}else{kk = kk+1;}
                                }
                                weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)])/(COHendo[inx((kk),y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)]);  // here kk-1 > kk
                                if(weight < 0){printf("MISTAKE, weight < 0");getchar();}
                                V_J3[(k-i_low)] = inter1d(weight,Vendo[inx((kk - 1),y,a,z,t)],Vendo[inx((kk),y,a,z,t)]);
                            }
                            
                        } // case 1 end
                        
                        
                        // CASE WHERE IN REGION 2 //
                        if(k > (i-1) && k < (i_low_up)){
                            
                            // VALUE ON SEGMENT 1 IS INTERPOLATED //
                            if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i-1),y,a,z,t)]){
                                slope = (Vendo[inx((i-1),y,a,z,t)] - Vendo[inx((i-2),y,a,z,t)])/(COHendo[inx((i-1),y,a,z,t)] - COHendo[inx((i-2),y,a,z,t)]);
                                V_J1[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i-1),y,a,z,t)]) + Vendo[inx((i-1),y,a,z,t)];
                            }
                            if(COHendo[inx(k,y,a,z,t)] < COHendo[inx(i_low,y,a,z,t)]){
                                slope = (Vendo[inx((i_low+1),y,a,z,t)] - Vendo[inx(i_low,y,a,z,t)])/(COHendo[inx((i_low+1),y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)]);
                                V_J1[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)]) + Vendo[inx((i_low),y,a,z,t)];
                            }
                            if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_low),y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i-1),y,a,z,t)]){ // INTERPOLATED
                                kk = (i_low+1);
                                while(COHendo[inx(k,y,a,z,t)] > COHendo[inx(kk,y,a,z,t)]){
                                    if(kk == (0)){printf("mistake, impossible case 0145");getchar();}else{kk = kk+1;}
                                }
                                weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)])/(COHendo[inx((kk),y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)]);  // here kk-1 > kk
                                if(weight < 0){printf("MISTAKE, weight < 0");getchar();}
                                V_J1[(k-i_low)] = inter1d(weight,Vendo[inx((kk-1),y,a,z,t)],Vendo[inx((kk),y,a,z,t)]);
                            }
                            
                            
                            // VALUE ON SEGMENT 3 IS INTERPOLATED //
                            if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_up),y,a,z,t)]){
                                slope = (Vendo[inx(i_up,y,a,z,t)] - Vendo[inx((i_up-1),y,a,z,t)])/(COHendo[inx((i_up),y,a,z,t)] - COHendo[inx((i_up-1),y,a,z,t)]);
                                V_J3[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i_up),y,a,z,t)]) + Vendo[inx((i_up),y,a,z,t)];
                            }
                            if(COHendo[inx(k,y,a,z,t)] < COHendo[inx((i_low_up),y,a,z,t)]){
                                slope = (Vendo[inx((i_low_up+1),y,a,z,t)] - Vendo[inx((i_low_up),y,a,z,t)])/(COHendo[inx((i_low_up+1),y,a,z,t)] - COHendo[inx((i_low_up),y,a,z,t)]);
                                V_J3[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i_low_up),y,a,z,t)]) + Vendo[inx((i_low_up),y,a,z,t)];
                            }
                            if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_low_up),y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i_up),y,a,z,t)]){ // INTERPOLATED
                                kk = (i_low_up+1);
                                while(COHendo[inx(k,y,a,z,t)] > COHendo[inx(kk,y,a,z,t)]){
                                    if(kk == (0)){printf("mistake, impossible case 0145");getchar();}else{kk = kk+1;}
                                }
                                weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)])/(COHendo[inx((kk),y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)]);  // here kk-1 > kk
                                if(weight < 0){printf("MISTAKE, weight < 0");getchar();}
                                V_J3[(k-i_low)] = inter1d(weight,Vendo[inx((kk-1),y,a,z,t)],Vendo[inx((kk),y,a,z,t)]);
                            }
                            
                        } // case 2 end
                        
                        
                        
                        // CASE WHERE IN REGION 3 //
                        if(k >= (i_low_up)){
                            
                            // VALUE ON SEGMENT 1 IS INTERPOLATED //
                            if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i-1),y,a,z,t)]){
                                slope = (Vendo[inx((i-1),y,a,z,t)] - Vendo[inx((i-2),y,a,z,t)])/(COHendo[inx((i-1),y,a,z,t)] - COHendo[inx((i-2),y,a,z,t)]);
                                V_J1[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i-1),y,a,z,t)]) + Vendo[inx((i-1),y,a,z,t)];
                            }
                            if(COHendo[inx(k,y,a,z,t)] < COHendo[inx(i_low,y,a,z,t)]){
                                slope = (Vendo[inx((i_low+1),y,a,z,t)] - Vendo[inx(i_low,y,a,z,t)])/(COHendo[inx((i_low+1),y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)]);
                                V_J1[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)]) + Vendo[inx((i_low),y,a,z,t)];
                            }
                            if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_low),y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i-1),y,a,z,t)]){ // INTERPOLATED
                                kk = (i_low+1);
                                while(COHendo[inx(k,y,a,z,t)] > COHendo[inx(kk,y,a,z,t)]){
                                    if(kk == (0)){printf("mistake, impossible case 0145");getchar();}else{kk = kk+1;}
                                }
                                weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)])/(COHendo[inx((kk),y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)]);  // here kk-1 > kk
                                if(weight < 0){printf("MISTAKE, weight < 0");getchar();}
                                V_J1[(k-i_low)] = inter1d(weight,Vendo[inx((kk-1),y,a,z,t)],Vendo[inx((kk),y,a,z,t)]);
                            }
                            
                            
                            // VALUE ON SEGMENT 3 IS INTERPOLATED //
                            V_J3[(k-i_low)] = Vendo[inx(k,y,a,z,t)];

                            
                        } // case 3 end
                        
                    
                    } // end for loop
                    
 
                
//                    /** CHECK **/
//                    printf("BOUND of M, low = %d, low_up = %d and up = %d", i_low, i_low_up, i_up);getchar();
//                    printf("CHECK COH: %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n", (i-5), (i-4), (i-3), (i-2), (i-1), (i), (i+1), (i+2), (i+3), (i+4), (i+5), (i+6), (i+7), (i+8));
//                    printf("CHECK COH: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", COHendo[inx((i-5),y,a,z,t)], COHendo[inx((i-4),y,a,z,t)],COHendo[inx((i-3),y,a,z,t)],COHendo[inx((i-2),y,a,z,t)],COHendo[inx((i-1),y,a,z,t)],COHendo[inx((i),y,a,z,t)],COHendo[inx((i+1),y,a,z,t)],COHendo[inx((i+2),y,a,z,t)],COHendo[inx((i+3),y,a,z,t)],COHendo[inx((i+4),y,a,z,t)],COHendo[inx((i+5),y,a,z,t)],COHendo[inx((i+6),y,a,z,t)],COHendo[inx((i+7),y,a,z,t)],COHendo[inx((i+8),y,a,z,t)]);
//                    printf("CHECK V: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", Vendo[inx((i-5),y,a,z,t)], Vendo[inx((i-4),y,a,z,t)],Vendo[inx((i-3),y,a,z,t)],Vendo[inx((i-2),y,a,z,t)],Vendo[inx((i-1),y,a,z,t)],Vendo[inx((i),y,a,z,t)],Vendo[inx((i+1),y,a,z,t)],Vendo[inx((i+2),y,a,z,t)],Vendo[inx((i+3),y,a,z,t)],Vendo[inx((i+4),y,a,z,t)],Vendo[inx((i+5),y,a,z,t)],Vendo[inx((i+6),y,a,z,t)],Vendo[inx((i+7),y,a,z,t)],Vendo[inx((i+8),y,a,z,t)]); getchar();
                    

//                    for(k = (i_low); k<(i_up+1); k++){
//                        if(V_J3[(k-i_low)] < -10000){
//                            printf("%d %d %d %d %d", i_low, i_up, i, ii, i_low_up); getchar();
//                    printf("BOUND of M = %d and %d", i_low, i_up);getchar();
//                    printf("CHECK COH: %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n", (i-5), (i-4), (i-3), (i-2), (i-1), (i), (i+1), (i+2), (i+3), (i+4), (i+5), (i+6), (i+7), (i+8));
//                    printf("CHECK COH: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", COHendo[inx((i-5),y,a,z,t)], COHendo[inx((i-4),y,a,z,t)],COHendo[inx((i-3),y,a,z,t)],COHendo[inx((i-2),y,a,z,t)],COHendo[inx((i-1),y,a,z,t)],COHendo[inx((i),y,a,z,t)],COHendo[inx((i+1),y,a,z,t)],COHendo[inx((i+2),y,a,z,t)],COHendo[inx((i+3),y,a,z,t)],COHendo[inx((i+4),y,a,z,t)],COHendo[inx((i+5),y,a,z,t)],COHendo[inx((i+6),y,a,z,t)],COHendo[inx((i+7),y,a,z,t)],COHendo[inx((i+8),y,a,z,t)]);
//                    printf("CHECK V: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", Vendo[inx((i-5),y,a,z,t)], Vendo[inx((i-4),y,a,z,t)],Vendo[inx((i-3),y,a,z,t)],Vendo[inx((i-2),y,a,z,t)],Vendo[inx((i-1),y,a,z,t)],Vendo[inx((i),y,a,z,t)],Vendo[inx((i+1),y,a,z,t)],Vendo[inx((i+2),y,a,z,t)],Vendo[inx((i+3),y,a,z,t)],Vendo[inx((i+4),y,a,z,t)],Vendo[inx((i+5),y,a,z,t)],Vendo[inx((i+6),y,a,z,t)],Vendo[inx((i+7),y,a,z,t)],Vendo[inx((i+8),y,a,z,t)]); getchar();
//
//                        }
//
//                        printf("k= %d %d, V1: %f, V3: %f \n", k,(k-i_low), V_J1[(k-i_low)], V_J3[(k-i_low)]);
//                    }
                    
                    
                    /** 5. DELETE SUBOPTIMAL POINTS **/
                    for(k = (i_low); k<(i_up+1); k++){
                        if(k >= i_low && k <= (i-1)){
                            if(V_J1[(k-i_low)] < V_J3[(k-i_low)]){COHendo[inx(k,y,a,z,t)] = -10000;}
                        }
                        if(k > (i-1) && k < (i_low_up)){
                            V_temp = max(V_J3[(k-i_low)], V_J1[(k-i_low)]);
                            if(Vendo[inx(k,y,a,z,t)] < V_temp){COHendo[inx(k,y,a,z,t)] = -10000;}
                            if(Vendo[inx(k,y,a,z,t)] > V_temp){printf("aie aie aie");getchar();}
                        }
                        if(k >= (i_low_up)){
                            if(V_J3[(k-i_low)] < V_J1[(k-i_low)]){COHendo[inx(k,y,a,z,t)] = -10000;}
                        }
                    }
                    
//                    for(k = (i_low); k<(i_up+1); k++){
//                        if(V_J3[(k-i_low)] < -10000){
//                            printf("%d %d %d %d %d", i_low, i_up, i, ii, i_low_up); getchar();
//                    printf("BOUND of M = %d and %d", i_low, i_up);getchar();
//                    printf("CHECK COH: %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n", (i-5), (i-4), (i-3), (i-2), (i-1), (i), (i+1), (i+2), (i+3), (i+4), (i+5), (i+6), (i+7), (i+8));
//                    printf("CHECK COH: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", COHendo[inx((i-5),y,a,z,t)], COHendo[inx((i-4),y,a,z,t)],COHendo[inx((i-3),y,a,z,t)],COHendo[inx((i-2),y,a,z,t)],COHendo[inx((i-1),y,a,z,t)],COHendo[inx((i),y,a,z,t)],COHendo[inx((i+1),y,a,z,t)],COHendo[inx((i+2),y,a,z,t)],COHendo[inx((i+3),y,a,z,t)],COHendo[inx((i+4),y,a,z,t)],COHendo[inx((i+5),y,a,z,t)],COHendo[inx((i+6),y,a,z,t)],COHendo[inx((i+7),y,a,z,t)],COHendo[inx((i+8),y,a,z,t)]);
//                    printf("CHECK V: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", Vendo[inx((i-5),y,a,z,t)], Vendo[inx((i-4),y,a,z,t)],Vendo[inx((i-3),y,a,z,t)],Vendo[inx((i-2),y,a,z,t)],Vendo[inx((i-1),y,a,z,t)],Vendo[inx((i),y,a,z,t)],Vendo[inx((i+1),y,a,z,t)],Vendo[inx((i+2),y,a,z,t)],Vendo[inx((i+3),y,a,z,t)],Vendo[inx((i+4),y,a,z,t)],Vendo[inx((i+5),y,a,z,t)],Vendo[inx((i+6),y,a,z,t)],Vendo[inx((i+7),y,a,z,t)],Vendo[inx((i+8),y,a,z,t)]); getchar();
//                    
//                        }
//                        
//                        printf("k= %d %d, V1: %f, V3: %f \n", k,(k-i_low), V_J1[(k-i_low)], V_J3[(k-i_low)]);
//                    }
                    
                
//                    /** CHECK **/
//                    printf("BOUND of M = %d and %d", i_low, i_up);getchar();
//                    printf("CHECK COH: %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n", (i-5), (i-4), (i-3), (i-2), (i-1), (i), (i+1), (i+2), (i+3), (i+4), (i+5), (i+6), (i+7), (i+8));
//                    printf("CHECK COH: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", COHendo[inx((i-5),y,a,z,t)], COHendo[inx((i-4),y,a,z,t)],COHendo[inx((i-3),y,a,z,t)],COHendo[inx((i-2),y,a,z,t)],COHendo[inx((i-1),y,a,z,t)],COHendo[inx((i),y,a,z,t)],COHendo[inx((i+1),y,a,z,t)],COHendo[inx((i+2),y,a,z,t)],COHendo[inx((i+3),y,a,z,t)],COHendo[inx((i+4),y,a,z,t)],COHendo[inx((i+5),y,a,z,t)],COHendo[inx((i+6),y,a,z,t)],COHendo[inx((i+7),y,a,z,t)],COHendo[inx((i+8),y,a,z,t)]);
//                    printf("CHECK V: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", Vendo[inx((i-5),y,a,z,t)], Vendo[inx((i-4),y,a,z,t)],Vendo[inx((i-3),y,a,z,t)],Vendo[inx((i-2),y,a,z,t)],Vendo[inx((i-1),y,a,z,t)],Vendo[inx((i),y,a,z,t)],Vendo[inx((i+1),y,a,z,t)],Vendo[inx((i+2),y,a,z,t)],Vendo[inx((i+3),y,a,z,t)],Vendo[inx((i+4),y,a,z,t)],Vendo[inx((i+5),y,a,z,t)],Vendo[inx((i+6),y,a,z,t)],Vendo[inx((i+7),y,a,z,t)],Vendo[inx((i+8),y,a,z,t)]); getchar();
                    
                    
                    /** 6. CHECK IF NO PROBLEM **/
                    for(k = (i_low+1); k<=(ii); k++){
                        coh_temp = COHendo[inx(k,y,a,z,t)];
                        coh_temp_next = COHendo[inx((k+1),y,a,z,t)];
                        
                        if(coh_temp_next != -10000 && coh_temp_next < coh_temp){printf("ERROR DC-EGM:: %d %f %f",k,coh_temp,coh_temp_next);getchar();}
                    }
                
                
                
                    /** 7. MARK DOWN IF KINK IS PRESENT **/
                    for(k = (i_low); k<=(i_up); k++){
                        coh_temp = COHendo[inx(k,y,a,z,t)];
                        coh_temp_next = COHendo[inx((k+1),y,a,z,t)];
                        
                        if(coh_temp == -10000 && coh_temp_next > coh_temp){kink_point[inx((k+1),y,a,z,t)] = 1;}
                        if(coh_temp > coh_temp_next && coh_temp_next == -10000){kink_point[inx((k),y,a,z,t)] = 2;}
                        
                        // here possibly we can have the case where 2 becomes 1 //
                    }
                
//                    for(k = i_low; k<(i_up+1); k++){
//                        printf("k= %d %d, TRUE V: %f,  V1: %f, V3: %f , kink = %d \n", k,(k-i_low), Vendo[inx(k,y,a,z,t)], V_J1[(k-i_low)], V_J3[(k-i_low)], kink_point[inx((k+1),y,a,z,t)]);
//                    }
                    
                } // end condition i_low_up != ii
                    
                //if(a == 43){printf("HERE2");getchar();}

                
                if(i_low_up == ii){ // case where there is only one kink region
                

                /** 4. Interpolate COHendo on the different path of V **/
                /* interpolate points between (i-1) and (i_low) to values in (ii) and (i_up) and (i) to (ii) */
                double V_J1[(i_up-i_low+1)], V_J2[(i_up-i_low+1)], V_J3[(i_up-i_low+1)];
                // all the points have to be computed here //

                for(k = (i_low); k<(i_up+1); k++){
                
                    // K iS IN SEGMENT 1 //
                    if(k < i && k >= i_low){
                    
                        /* ------ KNOW TRUE VALUE IN SEGMENT J1 ------ */
                        V_J1[(k-i_low)] = Vendo[inx(k,y,a,z,t)];    // true value
                        
                        
                        /* ------ INTERPOLATE K IN SEGMENT J2 ------ */
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i-1),y,a,z,t)]){ // case where you need to extrapolate above
                            printf("mistake, impossible case 02"); getchar(); // shouldn't happen
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] < COHendo[inx((ii),y,a,z,t)]){
                            slope = (Vendo[inx((ii-1),y,a,z,t)] - Vendo[inx((ii),y,a,z,t)])/(COHendo[inx((ii-1),y,a,z,t)] - COHendo[inx((ii),y,a,z,t)]);
                            V_J2[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((ii),y,a,z,t)]) + Vendo[inx((ii),y,a,z,t)];
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx(ii,y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i-1),y,a,z,t)]){   // normal interpolation in that case
                            kk = (i);
                            while(COHendo[inx(k,y,a,z,t)] < COHendo[inx(kk,y,a,z,t)]){
                                if(kk == (ii)){printf("mistake, impossible case 01");getchar();}else{kk = kk+1;}
                            }
                            weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk),y,a,z,t)])/(COHendo[inx((kk-1),y,a,z,t)] - COHendo[inx((kk),y,a,z,t)]);  // here kk-1 > kk
                            if(weight < 0){printf("MISTAKE, weight < 0");getchar();}
                            V_J2[(k-i_low)] = inter1d(weight,Vendo[inx((kk),y,a,z,t)],Vendo[inx((kk-1),y,a,z,t)]); // VERIF HERE
                            
                            //printf("w = %f, val1 = %f, val2 = %f, val = %f", weight, Vendo[inx((kk),y,a,z,t)],Vendo[inx((kk-1),y,a,z,t)], V_J2[(k-i_low)]); getchar();
                        }
                        
                        
                        /* ------ INTERPOLATE K IN SEGMENT J3 ------ */
                        if(COHendo[inx(k,y,a,z,t)] < COHendo[inx(ii,y,a,z,t)]){ // case where you need to extrapolate
                            slope = (Vendo[inx((ii+1),y,a,z,t)] - Vendo[inx((ii),y,a,z,t)])/(COHendo[inx((ii+1),y,a,z,t)] - COHendo[inx((ii),y,a,z,t)]);
                            V_J3[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx(ii,y,a,z,t)]) + Vendo[inx((ii),y,a,z,t)];
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_up),y,a,z,t)]){
                            printf("mistake, impossible case"); getchar();
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx(ii,y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i_up),y,a,z,t)]){   // normal interpolation in that case
                            kk = ii+1;
                            while(COHendo[inx(k,y,a,z,t)] > COHendo[inx(kk,y,a,z,t)]){
                                if(kk == i_up){printf("mistake, impossible case 2");getchar();}else{kk = kk+1;}
                            }
                            weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)])/(COHendo[inx(kk,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)]);
                            V_J3[(k-i_low)] = inter1d(weight,Vendo[inx((kk-1),y,a,z,t)],Vendo[inx((kk),y,a,z,t)]);
                        }
                    }   // END CONDITION OVER k IN SEGMENT 1 //
                    
                    
                    // K iS IN SEGMENT 2 //
                    if(k < (ii+1) && k >= (i-1)){
                        
                        /* ------ INTERPOLATE K IN SEGMENT J1 ------ */
                        if(COHendo[inx(k,y,a,z,t)] < COHendo[inx((i_low),y,a,z,t)]){ // can be possible at corner left
                            slope = (COHendo[inx((i_low+1),y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)])/(COHendo[inx((i_low+1),y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)]);
                            V_J1[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)]) + Vendo[inx((i_low),y,a,z,t)];

                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i-1),y,a,z,t)]){ // extrapolate
                            slope = (Vendo[inx((i-1),y,a,z,t)] - Vendo[inx((i-2),y,a,z,t)])/(COHendo[inx((i-1),y,a,z,t)] - COHendo[inx((i-2),y,a,z,t)]);
                            V_J1[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i-1),y,a,z,t)]) + Vendo[inx((i-1),y,a,z,t)];
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_low),y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i-1),y,a,z,t)]){   // normal interpolation in that case
                            kk = i_low+1;
                            while(COHendo[inx(k,y,a,z,t)] > COHendo[inx(kk,y,a,z,t)]){
                                if(kk == (i-1)){printf("mistake, impossible case 234");getchar();}else{kk = kk+1;}
                            }
                            weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)])/(COHendo[inx(kk,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)]);
                            V_J1[(k-i_low)] = inter1d(weight,Vendo[inx((kk-1),y,a,z,t)],Vendo[inx((kk),y,a,z,t)]);
                        }
                        
                        
                        /* ------ KNOW TRUE VALUE IN SEGMENT J2 ------ */
                        V_J2[(k-i_low)] = Vendo[inx(k,y,a,z,t)];    // true value
//                        if(k == 21){printf("%d %d %f %f %f %f",k,(k-i_low),  V_J2[(k-i_low)], Vendo[inx(k,y,a,z,t)], Vendo[inx((k-1),y,a,z,t)], Vendo[inx((k+1),y,a,z,t)]);getchar();}
                        
                        /* ------ INTERPOLATE K IN SEGMENT J3 ------ */
                        if(COHendo[inx(k,y,a,z,t)] < COHendo[inx(ii,y,a,z,t)]){ // case where you need to extrapolate
                            slope = (Vendo[inx((ii+1),y,a,z,t)] - Vendo[inx((ii),y,a,z,t)])/(COHendo[inx((ii+1),y,a,z,t)] - COHendo[inx((ii),y,a,z,t)]);
                            V_J3[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx(ii,y,a,z,t)]) + Vendo[inx((ii),y,a,z,t)];
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_up),y,a,z,t)]){
                            printf("mistake, impossible case 987"); getchar();
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx(ii,y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i_up),y,a,z,t)]){   // normal interpolation in that case
                            kk = ii+1;
                            while(COHendo[inx(k,y,a,z,t)] > COHendo[inx(kk,y,a,z,t)]){
                                if(kk == i_up){printf("mistake, impossible case 2");getchar();}else{kk = kk+1;}
                            }
                            weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)])/(COHendo[inx(kk,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)]);
                            V_J3[(k-i_low)] = inter1d(weight,Vendo[inx((kk-1),y,a,z,t)],Vendo[inx((kk),y,a,z,t)]);
                        }
                    }
                    
                    
                    // K iS IN SEGMENT 3 //
                    if(k < (i_up+1) && k >= ii){
                    
                        /* ------ INTERPOLATE K IN SEGMENT J1 ------ */
                        if(COHendo[inx(k,y,a,z,t)] < COHendo[inx((i_low),y,a,z,t)]){ // can be possible at left corner
                            slope = (COHendo[inx((i_low+1),y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)])/(COHendo[inx((i_low+1),y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)]);
                            V_J1[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i_low),y,a,z,t)]) + Vendo[inx((i_low),y,a,z,t)];
                            verifcase = 1;
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i-1),y,a,z,t)]){ // extrapolate
                            slope = (Vendo[inx((i-1),y,a,z,t)] - Vendo[inx((i-2),y,a,z,t)])/(COHendo[inx((i-1),y,a,z,t)] - COHendo[inx((i-2),y,a,z,t)]);
                            V_J1[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i-1),y,a,z,t)]) + Vendo[inx((i-1),y,a,z,t)];
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i_low),y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i-1),y,a,z,t)]){   // normal interpolation in that case
                            kk = i_low+1;
                            while(COHendo[inx(k,y,a,z,t)] > COHendo[inx(kk,y,a,z,t)]){
                                if(kk == (i-1)){printf("mistake, impossible case 234");getchar();}else{kk = kk+1;}
                            }
                            weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)])/(COHendo[inx(kk,y,a,z,t)] - COHendo[inx((kk-1),y,a,z,t)]);
                            V_J1[(k-i_low)] = inter1d(weight,Vendo[inx((kk-1),y,a,z,t)],Vendo[inx((kk),y,a,z,t)]);
                        }
                        
                        
                        /* ------ INTERPOLATE VALUE IN SEGMENT J2 ------ */
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx((i-1),y,a,z,t)]){ // case where you need to extrapolate above
                            slope = (Vendo[inx((i-1),y,a,z,t)] - Vendo[inx((i),y,a,z,t)])/(COHendo[inx((i-1),y,a,z,t)] - COHendo[inx((i),y,a,z,t)]); // (i-1) is greater than i
                            V_J2[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((i-1),y,a,z,t)]) + Vendo[inx((i-1),y,a,z,t)];
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] < COHendo[inx((ii),y,a,z,t)]){
                            slope = (Vendo[inx((ii-1),y,a,z,t)] - Vendo[inx((ii),y,a,z,t)])/(COHendo[inx((ii-1),y,a,z,t)] - COHendo[inx((ii),y,a,z,t)]);
                            V_J2[(k-i_low)] = slope*(COHendo[inx(k,y,a,z,t)] - COHendo[inx((ii),y,a,z,t)]) + Vendo[inx((ii),y,a,z,t)];
                        }
                        
                        if(COHendo[inx(k,y,a,z,t)] > COHendo[inx(ii,y,a,z,t)] && COHendo[inx(k,y,a,z,t)] < COHendo[inx((i-1),y,a,z,t)]){   // normal interpolation in that case
                            kk = (i);
                            while(COHendo[inx(k,y,a,z,t)] < COHendo[inx(kk,y,a,z,t)]){
                                if(kk == (ii)){printf("mistake, impossible case 03");getchar();}else{kk = kk+1;}
                            }
                            weight = (COHendo[inx(k,y,a,z,t)] - COHendo[inx((kk),y,a,z,t)])/(COHendo[inx((kk-1),y,a,z,t)] - COHendo[inx((kk),y,a,z,t)]);  // here kk-1 > kk
                            if(weight < 0){printf("MISTAKE, weight < 0");getchar();}
                            V_J2[(k-i_low)] = inter1d(weight,Vendo[inx((kk),y,a,z,t)],Vendo[inx((kk-1),y,a,z,t)]); // VERIF HERE
                            //printf("w = %f, val1 = %f, val2 = %f, val = %f", weight, Vendo[inx((kk),y,a,z,t)],Vendo[inx((kk-1),y,a,z,t)], V_J2[(k-i_low)]); getchar();
                        }
                        
                        
                        /* ------ KNOW TRUE VALUE OF K IN SEGMENT J3 ------ */
                        V_J3[(k-i_low)] = Vendo[inx(k,y,a,z,t)];    // true value
                        
                    }
                    
//                    if(k == 21){printf("%d %d %f %f %f %f",k,(k-i_low),  V_J2[(k-i_low)], Vendo[inx(k,y,a,z,t)], Vendo[inx((k-1),y,a,z,t)], Vendo[inx((k+1),y,a,z,t)]);getchar();}
                    
                }
                
//                printf("%f %f %f %f",  V_J2[(21-i_low)], Vendo[inx(21,y,a,z,t)], Vendo[inx((21-1),y,a,z,t)], Vendo[inx((21+1),y,a,z,t)]);getchar();




                
                

//                /** 5. FIND THE POINT OF THE KINK (using simple algebra when linear interpolation) **/
//                slope_J1 = (Vendo[inx((i-1),y,a,z,t)] - Vendo[inx((i-2),y,a,z,t)])/(COHendo[inx((i-1),y,a,z,t)] - COHendo[inx((i-2),y,a,z,t)]);
//                slope_J3 = (Vendo[inx((ii+1),y,a,z,t)] - Vendo[inx((ii),y,a,z,t)])/(COHendo[inx((ii+1),y,a,z,t)] - COHendo[inx((ii),y,a,z,t)]);
//                coh_kink = (-Vendo[inx(ii,y,a,z,t)] + Vendo[inx((i-1),y,a,z,t)] + COHendo[inx(ii,y,a,z,t)]*slope_J3 - COHendo[inx((i-1),y,a,z,t)]*slope_J1)/(-slope_J1 + slope_J3);
//                slope_J1 = (K[(i-1)] - Vendo[inx((i-2),y,a,z,t)])/(COHendo[inx((i-1),y,a,z,t)] - COHendo[inx((i-2),y,a,z,t)]);
//                slope_J3 = (Vendo[inx((ii+1),y,a,z,t)] - Vendo[inx((ii),y,a,z,t)])/(COHendo[inx((ii+1),y,a,z,t)] - COHendo[inx((ii),y,a,z,t)]);

                    
//                printf("COH KINKED: %f %f %f %f %f %f %f %f %f \n ", coh_kink, (coh_kink - COHendo[inx(ii,y,a,z,t)])*slope_J3 + Vendo[inx(ii,y,a,z,t)], (coh_kink - COHendo[inx((i-1),y,a,z,t)])*slope_J1 + Vendo[inx((i-1),y,a,z,t)], COHendo[inx(ii,y,a,z,t)], COHendo[inx((i-1),y,a,z,t)], Vendo[inx(ii,y,a,z,t)], Vendo[inx((i-1),y,a,z,t)], (coh_kink - COHendo[inx((i-1),y,a,z,t)])*slope_J1, (coh_kink - COHendo[inx(ii,y,a,z,t)])*slope_J3);
                
            
        
                
                /** 5. DELETE SUBOPTIMAL POINTS **/
                for(k = (i_low+1); k<(i_up); k++){
                    if(k < i && k >= i_low){
                        V_temp = max(V_J2[(k-i_low)], V_J3[(k-i_low)]);
                        if(V_J1[(k-i_low)] < V_temp){COHendo[inx(k,y,a,z,t)] = -10000;}
                    }
                    if(k < (ii+1) && k >= (i-1)){
                        V_temp = max(V_J3[(k-i_low)], V_J1[(k-i_low)]);
                        if(V_J2[(k-i_low)] < V_temp){COHendo[inx(k,y,a,z,t)] = -10000;}
                    }
                    if(k < (i_up+1) && k >= ii){
                        V_temp = max(V_J2[(k-i_low)], V_J1[(k-i_low)]);
                        if(V_J3[(k-i_low)] < V_temp){COHendo[inx(k,y,a,z,t)] = -10000;}
                    }
                }
                
                
                
//                * 7. REPLACE SUBOPTIMAL POINTS by KINKED VALUE *
                
                
                
                
                /** 6. CHECK IF NO PROBLEM **/
                for(k = (i_low+1); k<=(ii); k++){
                    coh_temp = COHendo[inx(k,y,a,z,t)];
                    coh_temp_next = COHendo[inx((k+1),y,a,z,t)];
                    
                    if(coh_temp_next != -10000 && coh_temp_next < coh_temp){printf("ERROR DC-EGM:: %d %f %f",k,coh_temp,coh_temp_next);getchar();}
                }
                

//                /** 7. MARK DOWN IF KINK IS PRESENT (IF NEXT STEP IS NOT PRESENT) **/
//                for(k = (i_low-1); k<=(i_up); k++){
//                    coh_temp = COHendo[inx(k,y,a,z,t)];
//                    coh_temp_next = COHendo[inx((k+1),y,a,z,t)];
//
//                    if(coh_temp == -10000 && coh_temp_next > coh_temp){kink_point[inx((k+1),y,a,z,t)] = 1;}
//                    if(coh_temp > coh_temp_next && coh_temp_next == -10000){kink_point[inx((k),y,a,z,t)] = 2;}
//
//                    // here possibly we can have the case where 2 becomes 1 //
//                }

                

                
//                if(verifcase == 1){
//                /** CHECK **/
//                printf("BOUND of M = %d and %d", i_low, i_up);getchar();
//
//                printf("CHECK COH: %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n", (i-5), (i-4), (i-3), (i-2), (i-1), (i), (i+1), (i+2), (i+3), (i+4), (i+5), (i+6), (i+7), (i+8));
//                printf("CHECK COH: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", COHendo[inx((i-5),y,a,z,t)], COHendo[inx((i-4),y,a,z,t)],COHendo[inx((i-3),y,a,z,t)],COHendo[inx((i-2),y,a,z,t)],COHendo[inx((i-1),y,a,z,t)],COHendo[inx((i),y,a,z,t)],COHendo[inx((i+1),y,a,z,t)],COHendo[inx((i+2),y,a,z,t)],COHendo[inx((i+3),y,a,z,t)],COHendo[inx((i+4),y,a,z,t)],COHendo[inx((i+5),y,a,z,t)],COHendo[inx((i+6),y,a,z,t)],COHendo[inx((i+7),y,a,z,t)],COHendo[inx((i+8),y,a,z,t)]);
//                printf("CHECK V: %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", Vendo[inx((i-5),y,a,z,t)], Vendo[inx((i-4),y,a,z,t)],Vendo[inx((i-3),y,a,z,t)],Vendo[inx((i-2),y,a,z,t)],Vendo[inx((i-1),y,a,z,t)],Vendo[inx((i),y,a,z,t)],Vendo[inx((i+1),y,a,z,t)],Vendo[inx((i+2),y,a,z,t)],Vendo[inx((i+3),y,a,z,t)],Vendo[inx((i+4),y,a,z,t)],Vendo[inx((i+5),y,a,z,t)],Vendo[inx((i+6),y,a,z,t)],Vendo[inx((i+7),y,a,z,t)],Vendo[inx((i+8),y,a,z,t)]); getchar();
//                }
//
//
//                if(verifcase == 1){
//                for(k = i_low; k<(i_up+1); k++){
//                    printf("k= %d %d, V1: %f, V2: %f, V3: %f , kink = %d \n", k,(k-i_low), V_J1[(k-i_low)], V_J2[(k-i_low)], V_J3[(k-i_low)], kink_point[inx((k),y,a,z,t)]);
//                }
//                }
                
                
                } // i_low_up == ii
                
                //if(a == 43){printf("HERE 3 %d %d %d %d", ii, i, i_low, i_up);getchar();}

                
                
                i = ii;
                
                
                //if(a == 43){printf("HERE last %d, %f %f",i, COHendo[inx(i,y,a,z,t)], COHendo[inx((i-1),y,a,z,t)]);getchar();}
                
            }else{ // END ICASE 1
                i = i + 1;
            }
            
            } // end i
        } // end z
    } // end t
} // end y


} // end UEC

