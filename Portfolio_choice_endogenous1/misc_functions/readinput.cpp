
// Function to read moment in a file // Transform a data file to a 1Darray // AG
void readinput(double *variable, const int dimension, const char *file) {    
    std::ifstream input(file); // Change the line to your file
    for (int i = 0; i < dimension; i++) {
        input >> variable[i];
    }
}

void readinputDOUBLE(double variable[], const int dimension, const char *file) {
    std::ifstream input(file); // Change the line to your file
    for (int i = 0; i < dimension; i++) {
        input >> variable[i];
    }
}


template<size_t dim2D_1, size_t dim2D_2>
void readINPUT2D(double (&variable)[dim2D_1][dim2D_2], const char *file) {
   
    FILE *temp_file;
    char scanval[80];

    temp_file = fopen(file, "r");
    for(int i = 0; i < dim2D_1; i++){
        for(int j = 0; j < dim2D_2; j++){
            fscanf(temp_file,"%s",scanval);
            variable[i][j] = atof(scanval);
        }
    }
    fclose(temp_file);

}



template<size_t dim1D>
void readINPUT1D(double (&variable)[dim1D], const char *file) {
    FILE *temp_file;
    char scanval[80];

    temp_file = fopen(file, "r");
    for(int i = 0; i < dim1D; i++){
        fscanf(temp_file,"%s",scanval);
        variable[i] = atof(scanval);
    }
    fclose(temp_file);
}
