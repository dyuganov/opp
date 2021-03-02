#include <random>
#include <iostream>
#include <time.h>

#include "Matrix.h"

// vectors duplicate in each process
void MpiV1NonlinearConjugateGradient(double* A, double* b, double* x, int rank, int size) {
    double startTime;
    if (rank == 0) startTime = MPI_Wtime();

    // algorithm variables
    auto* r = new double[N];
    auto* z = new double[N];

    double alpha = 0, beta = 0;
    auto* Az = new double[N];
    auto* Ax = new double[N];
    auto* alphaz = new double[N];
    auto* betaz = new double[N];
    auto* prev_r = new double[N];

    // MPI variables
    int matrixPartCapacity = N * N / size;
    int vectorPartCapacity = N / size;
    auto* matrixPart = new double[matrixPartCapacity];
    auto* mulResult = new double[vectorPartCapacity];

    if (rank == 0) {
        initRandMatrix(A);
        initRandVector(b);
        initRandVector(x);
    }

    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Scatter(A, matrixPartCapacity, MPI_DOUBLE, matrixPart, matrixPartCapacity, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    mulMatrixAndVector(matrixPart, N / size, x, mulResult);
    MPI_Allgather(mulResult, vectorPartCapacity, MPI_DOUBLE, Ax, vectorPartCapacity, MPI_DOUBLE, MPI_COMM_WORLD);

    subVector(b, Ax, r);
    memcpy(z, r, sizeof(double) * N);

    double prevDotProductR = dotProduct(r, r);
    double currDotProductR;
    const double bVectorLength = getVectorLength(b);
    double rVectorLength = getVectorLength(r);

    size_t iterationsCnt = 0;
    const double epsilon = 0.00001;

    // проверка на расхождение
    double prevEpsilonCheck = 0;
    size_t epsilonGrowCounter = 0;
    double epsilonCheck = rVectorLength / bVectorLength;

    while (epsilonCheck >= epsilon) {
        mulMatrixAndVector(matrixPart, N / size, z, mulResult);
        MPI_Allgather(mulResult, vectorPartCapacity, MPI_DOUBLE, Az, vectorPartCapacity, MPI_DOUBLE, MPI_COMM_WORLD);

        if(rank == 0){
            alpha = dotProduct(r, r) / dotProduct(Az, z);
            mulVectorScalar(z, alpha, alphaz);
            sumVector(x, alphaz, x);
            mulVectorScalar(Az, alpha, Az);
            memcpy(prev_r, r, sizeof(double) * N);
            subVector(r, Az, r);
            currDotProductR = dotProduct(r, r);
            beta = currDotProductR / prevDotProductR;
            prevDotProductR = currDotProductR;
            mulVectorScalar(z, beta, betaz);
            sumVector(r, betaz, z);

            //epsilonCheck = rVectorLength / bVectorLength;
            if(epsilonCheck > prevEpsilonCheck){
                ++epsilonGrowCounter;
                prevEpsilonCheck = epsilonCheck;
            }
            if(epsilonGrowCounter > 5){
                perror("Can't resolve the matrix.");
                std::cout << "Can't resolve the matrix." << std::endl;
                break;
            }

            ++iterationsCnt;
        }
        MPI_Bcast(r, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(z, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        rVectorLength = getVectorLength(r);
        epsilonCheck = rVectorLength / bVectorLength;
    }

    if(rank == 0) {
        std::cout << "Iterations: " << iterationsCnt << std::endl;
        std::cout << "TIME: " << MPI_Wtime() - startTime << std::endl;
    }

    delete[] alphaz;
    delete[] Az;
    delete[] Ax;
    delete[] betaz;
    delete[] prev_r;
    delete[] r;
    delete[] z;
}

int main(int argc, char* argv[]) {
    auto* A = new double[N * N];
    auto* b = new double[N];
    auto* x = new double[N];

    clock_t start = 0, end = 0;
    int size = 0, rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MpiV1NonlinearConjugateGradient(A, b, x, rank, size);
    MPI_Finalize();

    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
