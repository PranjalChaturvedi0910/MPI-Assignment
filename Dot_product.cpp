#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <vector>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int n = 1000000;
    int local_n = n / size;
    
    std::vector<double> local_a(local_n);
    std::vector<double> local_b(local_n);
    
    srand(time(NULL) + rank);
    for (int i = 0; i < local_n; i++) {
        local_a[i] = (double)rand() / RAND_MAX;
        local_b[i] = (double)rand() / RAND_MAX;
    }
    
    double local_dot = 0.0;
    for (int i = 0; i < local_n; i++) {
        local_dot += local_a[i] * local_b[i];
    }
    
    double global_dot;
    MPI_Reduce(&local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "Dot product of two vectors of size " << n << ": " << global_dot << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}