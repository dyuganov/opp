#include <iostream>
#include <windows.h>
//#include <mpi.h> // for cluster
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h" // for local use

// first (n1 x n2), second (n2 x n3), result (n1 x n3)
#define N1 (16)
#define N2 (16)
#define N3 (16)

// 16 cores
#define X_CORES_SIZE (8)
#define Y_CORES_SIZE (2)

#define X 0
#define Y 1

void printMatrix(const double* matrix, const int& n1, const int& n2){
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n2; ++j){
            std::cout << matrix[i * n1 + j] << " ";
        }
        std::cout << std::endl;
    }
}

void fillMatrix(double* matrix, const int& n1, const int& n2){
    for(int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            //matrix[j + i * n1] = rand() % 10;
            matrix[j + i * n1] = (i == j);
        }
    }
}


// A (n1 x n2), B (n2 x n3), result (n1 x n3)
void mulMatrixMatrix(const double* A, const double* B, double* result, const int& n1, const int& n2, const int& n3){
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n3; ++j){
            for(int k = 0; k < n2; ++k){
                //
            }
        }
    }
}


int main1(int argc, char* argv[]) {
    int size = 0, rank = 0;
    MPI_Comm colComm, rowComm, oldComm, newComm;
    oldComm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(oldComm, &size);
    MPI_Comm_rank(oldComm, &rank);
    if(size != X_CORES_SIZE * Y_CORES_SIZE) {
        if(rank == 0) fprintf(stderr, "Wrong cores num");
        MPI_Finalize();
        return 0;
    }

    int dims[2] = {X_CORES_SIZE, Y_CORES_SIZE};
    int coords[2] = {0};
    int periods[2] = {0};
    int remainDims[2] = {0};
    int mapping[X_CORES_SIZE * Y_CORES_SIZE] = {0};

    double* A = nullptr;
    double* B = nullptr;
    double* AB = nullptr;
    if (rank == 0) {
        A = new double[N1 * N2];
        B = new double[N2 * N3];
        AB = new double[N1 * N3];
        fillMatrix(A, N1, N2);
        fillMatrix(B, N2, N3);
        /*std::cout << "Matrix A" << std::endl;
        printMatrix(A, N1, N2);
        std::cout << std::endl << "Matrix B" << std::endl;
        printMatrix(B, N2, N3);*/
    }
    double* A_part = new double[(N1 * N2) / Y_CORES_SIZE];
    double* B_part = new double[(N2 * N3) / X_CORES_SIZE];
    double* AB_part = new double[(N1 * N3) / (X_CORES_SIZE * Y_CORES_SIZE)];

    double timeStart = 0;
    if(rank == 0) timeStart = MPI_Wtime();

    // main

    double timeFinish = 0;
    if(rank == 0) {
        timeFinish = MPI_Wtime();
        std::cout << "Time: " << timeFinish - timeStart << std::endl;
    }

    delete[] A;
    delete[] B;
    delete[] AB;
    delete[] A_part;
    delete[] B_part;
    delete[] AB_part;

    MPI_Finalize();
    return 0;
}
