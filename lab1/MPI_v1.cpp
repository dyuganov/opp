#include <random>
#include <iostream>

#include "Matrix.h"

void MpiV1NonlinearConjugateGradient(double* A, double* b, double* x, int rank, int size) {
    double startTime;
    if (rank == 0) startTime = MPI_Wtime();

    // algorithm variables
    double* r = new double[N];
    double* z = new double[N];

    double alpha = 0, beta = 0;
    double* Az = new double[N];
    double* Ax = new double[N];
    double* alphaz = new double[N];
    double* betaz = new double[N];
    double* prev_r = new double[N];

    // MPI variables
    int matrixPartSize = N * N / size;
    int vectorPartSize = N / size;
    double* matrixPart = new double[matrixPartSize];
    double* mulResult = new double[vectorPartSize];

    if (rank == 0) {
        initRandMatrix(A);
        initRandVector(b);
        initRandVector(x);
    }

    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Scatter(A, matrixPartSize, MPI_DOUBLE, matrixPart, matrixPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    mulMatrixAndVector(matrixPart, N / size, x, mulResult);
    MPI_Allgather(mulResult, vectorPartSize, MPI_DOUBLE, Ax, vectorPartSize, MPI_DOUBLE, MPI_COMM_WORLD);
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
        MPI_Allgather(mulResult, vectorPartSize, MPI_DOUBLE, Az, vectorPartSize, MPI_DOUBLE, MPI_COMM_WORLD);

        if(rank == 0){
            alpha = dotProduct(r, r) / dotProduct(Az, z);
            mulVectorScalar(z, alpha, alphaz);
            sumVector(x, alphaz, x);
            mulVectorScalar(Az, alpha, Az);

            subVector(r, Az, r);

            currDotProductR = dotProduct(r, r);
            beta = currDotProductR / prevDotProductR;
            prevDotProductR = currDotProductR;
            mulVectorScalar(z, beta, betaz);
            sumVector(r, betaz, z);

            if(epsilonCheck > prevEpsilonCheck){
                ++epsilonGrowCounter;
                prevEpsilonCheck = epsilonCheck;
            }
            if(epsilonGrowCounter > 5){
                fprintf(stderr, "Can't resolve the matrix.\n");
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

    int size = 0, rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) std::cout << "Name: " << argv[0] << std::endl << "ProcNum: " << size <<std::endl;

    MpiV1NonlinearConjugateGradient(A, b, x, rank, size);
    MPI_Finalize();

    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
