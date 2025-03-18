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
    std::vector<int> local_prefix_sum(local_n);
    
    srand(time(NULL) + rank);
    for (int i = 0; i < local_n; i++) {
        local_data[i] = rand() % 10;
    }
    
    int local_sum = 0;
    for (int i = 0; i < local_n; i++) {
        local_sum += local_data[i];
        local_prefix_sum[i] = local_sum;
    }
    
    int* prefix_sums = new int[size];
    MPI_Allgather(&local_sum, 1, MPI_INT, prefix_sums, 1, MPI_INT, MPI_COMM_WORLD);
    
    int offset = 0;
    for (int i = 0; i < rank; i++) {
        offset += prefix_sums[i];
    }
    
    for (int i = 0; i < local_n; i++) {
        local_prefix_sum[i] += offset;
    }
    
    std::vector<int> global_prefix_sum;
    if (rank == 0) {
        global_prefix_sum.resize(n);
    }
    
    MPI_Gather(local_prefix_sum.data(), local_n, MPI_INT, 
               global_prefix_sum.data(), local_n, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "First 10 elements of prefix sum: ";
        for (int i = 0; i < 10 && i < n; i++) {
            std::cout << global_prefix_sum[i] << " ";
        }
        std::cout << std::endl;
        
        std::cout << "Last 10 elements of prefix sum: ";
        for (int i = n - 10; i < n; i++) {
            std::cout << global_prefix_sum[i] << " ";
        }
        std::cout << std::endl;
    }
    
    delete[] prefix_sums;
    
    MPI_Finalize();
    return 0;
}