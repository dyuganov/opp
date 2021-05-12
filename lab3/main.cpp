#include <iostream>

//#include <mpi.h> // for cluster
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h" // for local use

// first (n1 x n2), second (n2 x n3), result (n1 x n3)
#define N1 (1000)
#define N2 (1000)
#define N3 (1000)

// 16 cores
#define X_COMP_SIZE (8)
#define Y_COMP_SIZE (2)

#define X 0
#define Y 1

struct Coordinates{
    int x = 0;
    int y = 0;
};

void printMatrix(const double* matrix, const int& n1, const int& n2){
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n2; ++j){
            std::cout << matrix[i * n1 + j] << " ";
        }
        std::cout << std::endl;
    }
}

void initMatrix(double* matrix, const int& n1, const int& n2){
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n2; ++j){
            matrix[j + i * n1] = rand()%10;
        }
    }
}

// A (n1 x n2), B (n2 x n3), result (n1 x n3)
void mulMatrixMatrix(const double* A, const double* B, double* result, const int& n1, const int& n2, const int& n3){
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n3; ++j){
            for(int k = 0; k < n2; ++k){

            }
        }
    }
}

int main(int argc, char* argv[]) {
    int size = 0, rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(size != X_COMP_SIZE * Y_COMP_SIZE) {
        fprintf(stderr, "Wrong cores num");
        MPI_Finalize();
        return 0;
    }

    



    MPI_Finalize();
    return 0;
}
