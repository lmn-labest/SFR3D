#include<ParallelMpi.h>

void mpiInit(int *argc,char **argv){

#ifdef _MPICH_
  MPI_Init(argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, comm);
  MPI_Comm_size(*comm, nPrcs);
  MPI_Comm_rank(*comm, myId);
/*... sem mpi*/
#else
  *comm = 0;
  *nPrcs= 1;
  *myId = 0;
#endif 
/*...................................................................*/
}
