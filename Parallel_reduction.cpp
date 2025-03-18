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
    
    std::vector<int> local_data(local_n);
    
    srand(time(NULL) + rank);
    for (int i = 0; i < local_n; i++) {
        local_data[i] = rand() % 100;
    }
    
    int local_sum = 0;
    for (int i = 0; i < local_n; i++) {
        local_sum += local_data[i];
    }
    
    int global_sum;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "Sum of " << n << " elements: " << global_sum << std::endl;
    }
    
    int local_max = local_data[0];
    for (int i = 1; i < local_n; i++) {
        if (local_data[i] > local_max) {
            local_max = local_data[i];
        }
    }
    
    int global_max;
    MPI_Reduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "Maximum value: " << global_max << std::endl;
    }
    
    int local_min = local_data[0];
    for (int i = 1; i < local_n; i++) {
        if (local_data[i] < local_min) {
            local_min = local_data[i];
        }
    }
    
    int global_min;
    MPI_Reduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "Minimum value: " << global_min << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}