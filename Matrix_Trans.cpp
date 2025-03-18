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
    
    int n = 1000;
    int local_rows = n / size;
    
    std::vector<int> local_matrix(local_rows * n);
    std::vector<int> local_transposed(n * local_rows);
    
    srand(time(NULL) + rank);
    for (int i = 0; i < local_rows; i++) {
        for (int j = 0; j < n; j++) {
            local_matrix[i * n + j] = rand() % 100;
        }
    }
    
    std::vector<int> sendcounts(size);
    std::vector<int> senddispls(size);
    std::vector<int> recvcounts(size);
    std::vector<int> recvdispls(size);
    
    for (int i = 0; i < size; i++) {
        sendcounts[i] = local_rows * (n / size);
        senddispls[i] = local_rows * (i * (n / size));
        recvcounts[i] = (n / size) * local_rows;
        recvdispls[i] = i * local_rows * (n / size);
    }
    
    std::vector<int> temp_matrix(local_rows * n);
    
    MPI_Alltoallv(local_matrix.data(), sendcounts.data(), senddispls.data(), MPI_INT,
                  temp_matrix.data(), recvcounts.data(), recvdispls.data(), MPI_INT,
                  MPI_COMM_WORLD);
    
    int block_size = n / size;
    for (int p = 0; p < size; p++) {
        for (int i = 0; i < block_size; i++) {
            for (int j = 0; j < local_rows; j++) {
                local_transposed[p * local_rows + j + i * n] = 
                    temp_matrix[p * block_size * local_rows + j * block_size + i];
            }
        }
    }
    
    if (rank == 0) {
        std::cout << "Matrix transposition complete." << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}