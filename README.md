# MPI Assignment

## Overview

This repository contains solutions for **MPI Assignments** which include parallel computing tasks using the **Message Passing Interface (MPI)**. The assignments cover fundamental MPI concepts and include various parallel algorithms.

## Contents

### **Assignment 2**
This part contains the following MPI-based programs:

1. **Monte Carlo Method for Estimating Pi** (`Monte_carlo.cpp`)
2. **Matrix Multiplication** (`Matrix_Mul.cpp`)
3. **Parallel Sorting using Odd-Even Sort** (`Odd_Even_sort.cpp`)
4. **Heat Distribution Simulation** (`Heat_Dist.cpp`)
5. **Parallel Reduction** (`Parallel_reduction.cpp`)
6. **Parallel Dot Product** (`Dot_product.cpp`)
7. **Parallel Prefix Sum (Scan)** (`prefix_sum.cpp`)
8. **Parallel Matrix Transposition** (`Matrix_Trans.cpp`)

### **Assignment 3**
This part contains the following MPI-based programs:

1. **DAXPY Computation (`Daxpy.cpp`)**  
   - Performs the **DAXPY** operation: `Y = a * X + Y` in parallel using **MPI Scatter and Gather**.  
   - Compares the execution time of sequential vs. parallel implementations.

2. **Pi Calculation using Numerical Integration (`Pi_calculation.cpp`)**  
   - Approximates the value of **π (pi)** using numerical integration.  
   - Implements **MPI Reduce** to sum up the results from multiple processes.

3. **Prime Number Calculation (`Prime_number.cpp`)**  
   - Computes **prime numbers** within a given range using parallel processing.  
   - Each process handles a subset of the range and collectively determines the prime numbers.

## Requirements

To compile and run these programs, you need:

- **MPI Library** (OpenMPI or MPICH)
- **C++ Compiler** (e.g., `g++`, `mpic++`)

## Installation

If you haven't installed MPI, install it using:

```sh
# Ubuntu/Debian
sudo apt install mpich

# macOS (using Homebrew)
brew install mpich
```

## How to Run

Each program can be compiled and executed using the following commands:

### **Assignment 2 Programs**

1. **Compile the program:**
   ```sh
   mpic++ filename.cpp -o filename.out
   ```
2. **Run the program with multiple processes:**
   ```sh
   mpirun -np 4 ./filename.out
   ```
   *(Change `4` to the desired number of processes.)*

#### **Example: Monte Carlo Pi Estimation**
```sh
mpic++ Monte_carlo.cpp -o Monte_carlo.out
mpirun -np 4 ./Monte_carlo.out
```

### **Assignment 3 Programs**

#### **DAXPY Computation**
```sh
mpic++ Daxpy.cpp -o Daxpy.out
mpirun -np 4 ./Daxpy.out
```

#### **Pi Calculation**
```sh
mpic++ Pi_calculation.cpp -o Pi_calculation.out
mpirun -np 4 ./Pi_calculation.out
```

#### **Prime Number Calculation**
```sh
mpic++ Prime_number.cpp -o Prime_number.out
mpirun -np 4 ./Prime_number.out
```

## Author

- **Pranjal Chaturvedi (102203290)**

## License

This project is for educational purposes. Modify and use it as needed!
