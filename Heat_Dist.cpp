#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>

#define N 100
#define MAX_ITERATIONS 1000
#define TOLERANCE 1e-6

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int rows_per_proc = N / size;
    int remainder = N % size;
    
    int local_rows = rows_per_proc + (rank < remainder ? 1 : 0);
    int start_row = rank * rows_per_proc + std::min(rank, remainder);
    
    std::vector<double> local_grid(local_rows * N, 0.0);
    std::vector<double> new_local_grid(local_rows * N, 0.0);
    
    if (rank == 0) {
        for (int j = 0; j < N; j++) {
            local_grid[j] = 100.0;
        }
    }
    
    if (rank == size - 1) {
        for (int j = 0; j < N; j++) {
            local_grid[(local_rows - 1) * N + j] = 0.0;
        }
    }
    
    for (int i = 0; i < local_rows; i++) {
        local_grid[i * N] = 100.0;
        local_grid[i * N + N - 1] = 100.0;
    }
    
    std::vector<double> top_row(N, 0.0);
    std::vector<double> bottom_row(N, 0.0);
    
    double global_diff = 1.0;
    int iteration = 0;
    
    while (global_diff > TOLERANCE && iteration < MAX_ITERATIONS) {
        int top_rank = (rank == 0) ? MPI_PROC_NULL : rank - 1;
        int bottom_rank = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
        
        MPI_Sendrecv(local_grid.data(), N, MPI_DOUBLE, top_rank, 0,
                     bottom_row.data(), N, MPI_DOUBLE, bottom_rank, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        MPI_Sendrecv(local_grid.data() + (local_rows - 1) * N, N, MPI_DOUBLE, bottom_rank, 0,
                     top_row.data(), N, MPI_DOUBLE, top_rank, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        double local_diff = 0.0;
        
        for (int i = 0; i < local_rows; i++) {
            for (int j = 0; j < N; j++) {
                if (i == 0 && j == 0) continue;
                if (i == 0 && j == N - 1) continue;
                if (i == local_rows - 1 && j == 0) continue;
                if (i == local_rows - 1 && j == N - 1) continue;
                if (j == 0 || j == N - 1) continue;
                if ((rank == 0 && i == 0) || (rank == size - 1 && i == local_rows - 1)) continue;
                
                double up = (i == 0) ? top_row[j] : local_grid[(i - 1) * N + j];
                double down = (i == local_rows - 1) ? bottom_row[j] : local_grid[(i + 1) * N + j];
                double left = local_grid[i * N + j - 1];
                double right = local_grid[i * N + j + 1];
                
                new_local_grid[i * N + j] = 0.25 * (up + down + left + right);
                
                local_diff = std::max(local_diff, std::abs(new_local_grid[i * N + j] - local_grid[i * N + j]));
            }
        }
        
        local_grid.swap(new_local_grid);
        
        MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        iteration++;
    }
    
    if (rank == 0) {
        std::cout << "Converged after " << iteration << " iterations" << std::endl;
        std::cout << "Final maximum difference: " << global_diff << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}