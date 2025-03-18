#include <iostream>
#include <mpi.h>
#include <vector>
#include <cstdlib>
#include <ctime>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    const int N = 1 << 16;  // 2^16
    const double a = 2.5;   // Scalar value
    
    int local_size = N / size;
    int start_idx = rank * local_size;
    int end_idx = start_idx + local_size;
    
    std::vector<double> X(N), Y(N);
    std::vector<double> local_X(local_size), local_Y(local_size);
    
    if (rank == 0) {
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            X[i] = (double)rand() / RAND_MAX;
            Y[i] = (double)rand() / RAND_MAX;
        }
        
        double start_time = MPI_Wtime();
        
        for (int i = 0; i < N; i++) {
            X[i] = a * X[i] + Y[i];
        }
        
        double end_time = MPI_Wtime();
        std::cout << "Sequential execution time: " << (end_time - start_time) << " seconds" << std::endl;
        
        for (int i = 0; i < N; i++) {
            X[i] = (double)rand() / RAND_MAX;
            Y[i] = (double)rand() / RAND_MAX;
        }
    }
    
    double parallel_start_time = MPI_Wtime();
    
    MPI_Scatter(X.data(), local_size, MPI_DOUBLE, local_X.data(), local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(Y.data(), local_size, MPI_DOUBLE, local_Y.data(), local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for (int i = 0; i < local_size; i++) {
        local_X[i] = a * local_X[i] + local_Y[i];
    }
    
    MPI_Gather(local_X.data(), local_size, MPI_DOUBLE, X.data(), local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    double parallel_end_time = MPI_Wtime();
    
    if (rank == 0) {
        std::cout << "Parallel execution time: " << (parallel_end_time - parallel_start_time) << " seconds" << std::endl;
        double speedup = (end_time - start_time) / (parallel_end_time - parallel_start_time);
        std::cout << "Speedup: " << speedup << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}