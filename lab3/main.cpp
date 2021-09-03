#include "stdio.h"
#include <iostream>
#include <cstring>

#ifdef __unix__
#include <mpi.h>
#elif defined(_WIN32) || defined(WIN32)
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#endif

#define N1 2000
#define N2 2000
#define N3 2000
#define X_CORES_NUM 8
#define Y_CORES_NUM 2
#define X 0
#define Y 1

void mulMatrix(const double* A, const double* B, double* AB, int n1, int n2, int n3){
    for (int i = 0; i < n1; ++i){
        double* ba = AB + i * n3;
        for (int j = 0; j < n3; ++j) ba[j] = 0;
        for (int k = 0; k < n2; ++k){
            const double* b = B + k * n3;
            double a = A[i * n2 + k];
            for (int j = 0; j < n3; ++j){
                ba[j] += a * b[j];
            }
        }
    }
}

void printMatrix(double* matrix, int n1, int n2){
    std::cout << std::endl;
    for (size_t i = 0; i < n1 * n2; ++i){
        if (i % n2 == 0){
            std::cout << std::endl;
        }
        std::cout << matrix[i] << " ";
    }
    std::cout << std::endl;
}

void generateMatrix(double* matrix, int n1, int n2){
    for (size_t i = 0; i < n1; ++i){
        for (size_t j = 0; j < n2; ++j){
            matrix[i * n2 + j] = rand() % 10;
        }
    }
}

void copyMatrixFromPart(double* A, const double* B, int n1, int n2){
    for (int i = 0; i < n1 / Y_CORES_NUM; ++i){
        for (int j = 0; j < n2 / X_CORES_NUM; ++j){
            A[i * n2 + j] = B[i * n2 / X_CORES_NUM + j];
        }
    }
}


void copyMatrixToPart(double* A, const double* B, int n1, int n2) {
    for (int i = 0; i < n1; ++i){
        for (int j = 0; j < n2 / X_CORES_NUM; ++j){
            A[i * n2 / X_CORES_NUM + j] = B[i * n2 + j];
        }
    }
}

int main(int argc, char** argv){
    int size = 0, rank = 0;
    int dims[2] = {X_CORES_NUM, Y_CORES_NUM};
    int periods[2] = { 0, 0 };
    int cords[2] = { 0 ,0 };
    int remainDims[2] = {0, 0 };
    int mapping[X_CORES_NUM * Y_CORES_NUM] = {0};
    double timeStart, timeFinish;
    MPI_Comm oldComm, newComm, rowComm, colComm;
    oldComm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(oldComm, &size);
    MPI_Comm_rank(oldComm, &rank);

    if (size != X_CORES_NUM * Y_CORES_NUM){
        fprintf(stderr, "size != X_CORES_NUM * Y_CORES_NUM\n");
        MPI_Finalize();
        return 0;
    }

    double* A = nullptr;
    double* B = nullptr;
    double* AB = nullptr;
    double* A_matrixPart = new double[(N1 * N2) / Y_CORES_NUM]();
    double* B_matrixPart = new double[(N2 * N3) / X_CORES_NUM]();
    double* AB_matrixPart = new double[(N1 * N3) / X_CORES_NUM * Y_CORES_NUM]();

    if (rank == 0){
        A = new double[N1 * N2]();
        B = new double[N2 * N3]();
        AB = new double[N1 * N3]();
        generateMatrix(A, N1, N2);
        generateMatrix(B, N2, N3);
    }

    if (rank == 0) timeStart = MPI_Wtime();

    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(oldComm, 2, dims, periods, 1, &newComm);
    MPI_Comm_rank(newComm, &rank);
    MPI_Cart_coords(newComm, rank, 2, cords);

    const int TAG1 = 30;
    if (rank != 0) MPI_Send(cords, 2, MPI_INT, 0, TAG1, newComm);
    if (rank == 0){
        for (int i = 1; i < size; ++i){
            MPI_Recv(cords, 2, MPI_INT, i, TAG1, newComm, MPI_STATUS_IGNORE);
            mapping[i] = cords[X] * 10 + cords[Y];
        }
        cords[X] = 0;
        cords[Y] = 0;
    }
    
    remainDims[Y] = 0;
    remainDims[X] = 1;
    MPI_Cart_sub(newComm, remainDims, &rowComm);

    remainDims[Y] = 1;
    remainDims[X] = 0;
    MPI_Cart_sub(newComm, remainDims, &colComm);

    if (cords[X] == 0){
        MPI_Scatter(A, (N1*N2) / Y_CORES_NUM, MPI_DOUBLE, A_matrixPart, (N1 * N2) / Y_CORES_NUM, MPI_DOUBLE, 0, colComm);
    }

    MPI_Bcast(A_matrixPart, (N1 * N2) / Y_CORES_NUM, MPI_DOUBLE, 0, rowComm);

    MPI_Datatype COLUMN;
    MPI_Type_vector(N2, N3 / X_CORES_NUM, N3, MPI_DOUBLE, &COLUMN);
    MPI_Type_commit(&COLUMN);

    const int TAG2 = 10;
    if (rank == 0){
        copyMatrixToPart(B_matrixPart, B, N2, N3);
        for (int i = 1; i < X_CORES_NUM; i++){
            MPI_Send(B + i * (N3 / X_CORES_NUM), 1, COLUMN, i, TAG2, rowComm);
        }
    }
    if (cords[Y] == 0 && rank != 0){
        MPI_Recv(B_matrixPart, N2 * (N3 / X_CORES_NUM), MPI_DOUBLE, 0, TAG2, rowComm, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(B_matrixPart, (N2 * N3) / X_CORES_NUM, MPI_DOUBLE, 0, colComm);

    mulMatrix(A_matrixPart, B_matrixPart, AB_matrixPart, N1 / Y_CORES_NUM, N2, N3 / X_CORES_NUM);

    MPI_Datatype MINOR;
    MPI_Type_vector(N1 / Y_CORES_NUM, N3 / X_CORES_NUM, N3, MPI_DOUBLE, &MINOR);
    MPI_Type_commit(&MINOR);

    const int TAG3 = 20;
    if (rank != 0){
        MPI_Send(AB_matrixPart, (N1 * N3) / (X_CORES_NUM * Y_CORES_NUM), MPI_DOUBLE, 0, TAG3, newComm);
    }
    if (rank == 0){
        for (int i = 1; i < size; ++i){
            cords[X] = mapping[i] / 10;
            cords[Y] = mapping[i] - cords[X] * 10;
            MPI_Recv(AB + cords[X] * (N3 / X_CORES_NUM) + cords[Y] * (N3 * (N1 / Y_CORES_NUM)), 1, MINOR, i, TAG3, newComm, MPI_STATUS_IGNORE);
        }
        copyMatrixFromPart(AB, AB_matrixPart, N1, N3);
    }

    if (rank == 0) {
        timeFinish = MPI_Wtime();
        std::cout << "Time: " << (timeFinish - timeStart) << std::endl;
    }

    delete[] B;
    delete[] AB;
    delete[] A;
    delete[] A_matrixPart;
    delete[] B_matrixPart;
    delete[] AB_matrixPart;

    MPI_Type_free(&COLUMN);
    MPI_Type_free(&MINOR);
    MPI_Finalize();
    return 0;
}
