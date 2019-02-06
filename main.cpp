#include <iostream>
#include "Interacting_Harmonic_Oscillator.hh"

#include <omp.h>
#include <mpi/mpi.h>

int main(int argc, char * argv[]) {
    MPI_Init(&argc,&argv);
    int prank,psize;
    MPI_Comm_rank(MPI_COMM_WORLD,&prank);
    MPI_Comm_size(MPI_COMM_WORLD,&psize);

    std::cout << "Hello, World!" << std::endl;
   Interacting_Harmonic_Oscillator sim (64, 0.0025,
                                         10000, 300,
                                         0.01, 1,
                                         1000,
                                         100,
                                         "correlations"+to_string(8*prank)+".csv", 0.25,
                                         20,  1,  225,  1.486*10000,
                                          285, 903, 75,
                                          0.2, 5, 0.0005,
                                          75,  0.001,
                                          64, 8*prank);


   sim.sample();
   MPI_Finalize();

    return 0;
}