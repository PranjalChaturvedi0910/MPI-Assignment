#include <iostream>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <cmath>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    const int points_per_proc = 1000000;
    
    srand(time(NULL) + rank);
    
    int local_points_in_circle = 0;
    
    for (int i = 0; i < points_per_proc; i++) {
        double x = (double)rand() / RAND_MAX;
        double y = (double)rand() / RAND_MAX;
        
        if (x*x + y*y <= 1.0) {
            local_points_in_circle++;
        }
    }
    
    int total_points_in_circle = 0;
    MPI_Reduce(&local_points_in_circle, &total_points_in_circle, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        int total_points = points_per_proc * size;
        double pi_estimate = 4.0 * total_points_in_circle / total_points;
        
        std::cout << "Estimated Pi: " << pi_estimate << std::endl;
        std::cout << "Error: " << std::abs(pi_estimate - M_PI) << std::endl;
        std::cout << "Points used: " << total_points << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}