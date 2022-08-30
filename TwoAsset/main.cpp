#include <mpi.h>

int main(){
    MPI::Init(argc,argv);

    int num_proc = MPI::COMM_WORLD.Get_size();
    
}