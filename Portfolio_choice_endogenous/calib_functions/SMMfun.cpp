
// Function SMM Loss Function (GMM) // AG
double SMM(double *momentsobv, double *momentsgen, double covmatrixe[nb_moments][nb_moments]) {
	double *arg1, *matrix1, minLossfun;
	arg1 = (double *) calloc((nb_moments), sizeof(double));
	matrix1 = (double *) calloc((nb_moments), sizeof(double));
	
	for(int i=0; i < nb_moments; i++) {
		arg1[i] = fabs(momentsobv[i] - momentsgen[i])/momentsobv[i]; // measure the precision
        //printf("arg = %f | %f %f \n",arg1[i], momentsobv[i], momentsgen[i]);getchar();
	} // get a 1xN matrix
	
//    for(int i=0; i < nb_moments; i++) {
//        matrix1[i] = 0.0;
//        for(int j = 0; j < nb_moments; j++) {
//            matrix1[i] += arg1[j]*covmatrixe[j][i]; // covmatrix is of size NxN
//        }
//    } // get a 1xN matrix
	
	minLossfun = 0.0;
	for(int i=0; i < nb_moments; i++) {
		//minLossfun += matrix1[i]*arg1[i];
        minLossfun += arg1[i];
	} // get a 1x1 object
	
	return minLossfun;
}

