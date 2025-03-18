#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <omp.h>

#define N 70

void initialize_matrix(double* matrix) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i * N + j] = (double)rand() / RAND_MAX;
        }
    }
}

void sequential_multiply(double* A, double* B, double* C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i * N + j] = 0.0;
            for (int k = 0; k < N; k++) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    double* A = nullptr;
    double* B = nullptr;
    double* C = nullptr;
    double* C_seq = nullptr;

    if (rank == 0) {
        srand(time(NULL));
        A = new double[N * N];
        B = new double[N * N];
        C = new double[N * N];
        C_seq = new double[N * N];
        
        initialize_matrix(A);
        initialize_matrix(B);
        
        std::cout << "Matrix size: " << N << "x" << N << std::endl;
        
        double start_time = omp_get_wtime();
        sequential_multiply(A, B, C_seq);
        double run_time = omp_get_wtime() - start_time;
        
        std::cout << "Sequential time: " << run_time << " seconds" << std::endl;
    }
    
    double* B_local = new double[N * N];
    MPI_Bcast(B_local, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int rows_per_proc = N / size;
    int remainder = N % size;
    
    int local_rows = rows_per_proc + (rank < remainder ? 1 : 0);
    
    int* sendcounts = new int[size];
    int* displs = new int[size];
    
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (rows_per_proc + (i < remainder ? 1 : 0)) * N;
        displs[i] = (i > 0) ? displs[i-1] + sendcounts[i-1] : 0;
    }
    
    double* A_local = new double[local_rows * N];
    double* C_local = new double[local_rows * N];
    
    double start_time_parallel = 0.0;
    if (rank == 0) {
        start_time_parallel = omp_get_wtime();
    }
    
    MPI_Scatterv(A, sendcounts, displs, MPI_DOUBLE, A_local, local_rows * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for (int i = 0; i < local_rows; i++) {
        for (int j = 0; j < N; j++) {
            C_local[i * N + j] = 0.0;
            for (int k = 0; k < N; k++) {
                C_local[i * N + j] += A_local[i * N + k] * B_local[k * N + j];
            }
        }
    }
    
    MPI_Gatherv(C_local, local_rows * N, MPI_DOUBLE, C, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        double run_time_parallel = omp_get_wtime() - start_time_parallel;
        std::cout << "Parallel time: " << run_time_parallel << " seconds" << std::endl;
        std::cout << "Speedup: " << (run_time_parallel > 0 ? run_time_parallel / run_time_parallel : 0) << std::endl;
    }
    
    delete[] A_local;
    delete[] B_local;
    delete[] C_local;
    delete[] sendcounts;
    delete[] displs;
    
    if (rank == 0) {
        delete[] A;
        delete[] B;
        delete[] C;
        delete[] C_seq;
    }
    
    MPI_Finalize();
    return 0;
}