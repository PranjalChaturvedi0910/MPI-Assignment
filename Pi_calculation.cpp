#include <iostream>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    static long num_steps = 100000;
    double step;
    double pi, sum = 0.0;
    
    MPI_Bcast(&num_steps, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    
    step = 1.0 / (double)num_steps;
    
    int steps_per_process = num_steps / size;
    int start = rank * steps_per_process;
    int end = (rank == size - 1) ? num_steps : start + steps_per_process;
    
    double x;
    double local_sum = 0.0;
    
    for (int i = start; i < end; i++) {
        x = (i + 0.5) * step;
        local_sum += 4.0 / (1.0 + x * x);
    }
    
    MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        pi = step * sum;
        std::cout << "Pi approximation: " << pi << std::endl;
        std::cout << "Error: " << std::abs(pi - 3.14159265358979323846) << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}