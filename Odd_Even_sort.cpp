#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <algorithm>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int n = 1000;
    int local_n = n / size;
    
    int* global_data = nullptr;
    if (rank == 0) {
        global_data = new int[n];
        srand(time(NULL));
        for (int i = 0; i < n; i++) {
            global_data[i] = rand() % 10000;
        }
    }
    
    int* local_data = new int[local_n];
    
    MPI_Scatter(global_data, local_n, MPI_INT, local_data, local_n, MPI_INT, 0, MPI_COMM_WORLD);
    
    std::sort(local_data, local_data + local_n);
    
    int* recv_buf = new int[local_n];
    
    for (int i = 0; i < size; i++) {
        if (i % 2 == 0) {
            if (rank % 2 == 0 && rank < size - 1) {
                MPI_Send(local_data, local_n, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                MPI_Recv(recv_buf, local_n, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                merge_arrays(local_data, recv_buf, local_n);
            } else if (rank % 2 == 1) {
                MPI_Recv(recv_buf, local_n, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(local_data, local_n, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
                merge_arrays(recv_buf, local_data, local_n);
            }
        } else {
            if (rank % 2 == 1 && rank < size - 1) {
                MPI_Send(local_data, local_n, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                MPI_Recv(recv_buf, local_n, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                merge_arrays(local_data, recv_buf, local_n);
            } else if (rank % 2 == 0 && rank > 0) {
                MPI_Recv(recv_buf, local_n, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(local_data, local_n, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
                merge_arrays(recv_buf, local_data, local_n);
            }
        }
    }
    
    MPI_Gather(local_data, local_n, MPI_INT, global_data, local_n, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        bool is_sorted = true;
        for (int i = 0; i < n - 1; i++) {
            if (global_data[i] > global_data[i + 1]) {
                is_sorted = false;
                break;
            }
        }
        std::cout << "Array is " << (is_sorted ? "sorted" : "not sorted") << std::endl;
        delete[] global_data;
    }
    
    delete[] local_data;
    delete[] recv_buf;
    
    MPI_Finalize();
    return 0;
}

void merge_arrays(int* a, int* b, int n) {
    int* merged = new int[2 * n];
    int i = 0, j = 0, k = 0;
    
    while (i < n && j < n) {
        if (a[i] <= b[j]) {
            merged[k++] = a[i++];
        } else {
            merged[k++] = b[j++];
        }
    }
    
    while (i < n) {
        merged[k++] = a[i++];
    }
    
    while (j < n) {
        merged[k++] = b[j++];
    }
    
    for (i = 0; i < n; i++) {
        a[i] = merged[i];
    }
    
    delete[] merged;
}