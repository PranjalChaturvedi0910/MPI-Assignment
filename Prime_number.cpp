#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>

bool is_prime(int n) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    
    int i = 5;
    while (i * i <= n) {
        if (n % i == 0 || n % (i + 2) == 0) return false;
        i += 6;
    }
    return true;
}

void master_process(int size, int max_value) {
    std::vector<int> primes;
    int active_slaves = size - 1;
    int next_number = 2;
    int number_received;
    MPI_Status status;
    
    while (active_slaves > 0) {
        MPI_Recv(&number_received, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        
        int slave_rank = status.MPI_SOURCE;
        
        if (number_received == 0) {
            // Slave is just starting
        } else if (number_received > 0) {
            // Received a prime number
            primes.push_back(number_received);
        }
        
        if (next_number <= max_value) {
            MPI_Send(&next_number, 1, MPI_INT, slave_rank, 0, MPI_COMM_WORLD);
            next_number++;
        } else {
            // Send termination signal
            int terminate = -1;
            MPI_Send(&terminate, 1, MPI_INT, slave_rank, 0, MPI_COMM_WORLD);
            active_slaves--;
        }
    }
    
    std::cout << "Found " << primes.size() << " prime numbers up to " << max_value << std::endl;
    std::cout << "First 10 primes: ";
    for (int i = 0; i < 10 && i < primes.size(); i++) {
        std::cout << primes[i] << " ";
    }
    std::cout << std::endl;
    
    if (primes.size() > 10) {
        std::cout << "Last 10 primes: ";
        for (int i = primes.size() - 10; i < primes.size(); i++) {
            std::cout << primes[i] << " ";
        }
        std::cout << std::endl;
    }
}

void slave_process(int rank) {
    int number_to_test;
    int result;
    
    // Send initial request
    result = 0;
    MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    
    while (true) {
        MPI_Recv(&number_to_test, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        if (number_to_test == -1) {
            // Termination signal
            break;
        }
        
        if (is_prime(number_to_test)) {
            result = number_to_test;
        } else {
            result = -number_to_test;
        }
        
        MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int max_value = 10000;
    
    if (rank == 0) {
        master_process(size, max_value);
    } else {
        slave_process(rank);
    }
    
    MPI_Finalize();
    return 0;
}