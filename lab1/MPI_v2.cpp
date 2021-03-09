#include <random>
#include <iostream>
#include <time.h>

#include "Matrix.h"

// vectors duplicate in each process
void MpiV2NonlinearConjugateGradient(double* A, double* b, double* x, int rank, int size) {
    double startTime;
    if (rank == 0) startTime = MPI_Wtime();

    // algorithm variables
    auto* r = new double[N]; // на 0 потоке только?
    auto* z = new double[N];

    double alpha = 0, beta = 0;
    //auto* Az = new double[N]; // не выделять память?
    //auto* Ax = new double[N]; // сейчас она во всех потоках
    //auto* alphaz = new double[N];
    //auto* betaz = new double[N];
    //auto* prev_r = new double[N];

    // MPI variables
    int matrixPartCapacity = N * N / size;
    int vectorPartCapacity = N / size;
    auto* matrixPart = new double[matrixPartCapacity];
    auto* mulMatrixAndVectorResult = new double[vectorPartCapacity];
    //auto* vectorPart = new double[vectorPartCapacity];

    auto* Az_vecPart = new double[vectorPartCapacity];
    auto* alphaz_vecPart = new double[vectorPartCapacity];
    auto* betaz_vecPart = new double[vectorPartCapacity];
    auto* b_vecPart = new double[vectorPartCapacity];
    auto* x_vecPart = new double[vectorPartCapacity];
    auto* z_vecPart = new double[vectorPartCapacity];

    /*auto* alphaz_Part = new double[vectorPartCapacity];
    auto* AzPart = new double[vectorPartCapacity];
    auto* betazPart = new double[vectorPartCapacity];*/


    if (rank == 0) {
        initRandMatrix(A);
        initRandVector(b);
        initRandVector(x);
    }

    //MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(x, vectorPartCapacity, MPI_DOUBLE, x_vecPart, vectorPartCapacity, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, vectorPartCapacity, MPI_DOUBLE, b_vecPart, vectorPartCapacity, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(A, matrixPartCapacity, MPI_DOUBLE, matrixPart, matrixPartCapacity, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    mulMatrixAndVectorParts(matrixPart, x_vecPart, Ax, size, rank);
    //mulMatrixAndVector(matrixPart, N / size, x, mulMatrixAndVectorResult);
    //MPI_Allgather(mulMatrixAndVectorResult, vectorPartCapacity, MPI_DOUBLE, Ax, vectorPartCapacity, MPI_DOUBLE, MPI_COMM_WORLD);

    auto* Ax_vecPart = new double[vectorPartCapacity];
    //MPI_Scatter(Ax, vectorPartCapacity, MPI_DOUBLE, Ax_vecPart, vectorPartCapacity, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    auto* r_vecPart = new double[vectorPartCapacity];
    subVector(b_vecPart, Ax_vecPart, r_vecPart, vectorPartCapacity);
    memcpy(z_vecPart, r_vecPart, sizeof(double) * vectorPartCapacity);

    //double prevDotProductR = dotProduct(r, r);
    double prevDotProductRPart = dotProduct(r_vecPart, r_vecPart, vectorPartCapacity);
    double currDotProductRPart;

    double b_VecPartLength = getVectorPartSqr(b_vecPart, vectorPartCapacity);
    double b_VecLength = 0;
    MPI_Allreduce(&b_VecPartLength, &b_VecLength, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    b_VecLength = sqrt(b_VecLength);
    //double rVectorLength = getVectorLength(r);
    double r_VectorPartLength = getVectorPartSqr(r_vecPart, vectorPartCapacity);
    double r_VectorLength = 0;
    MPI_Allreduce(&r_VectorPartLength, &r_VectorLength, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    r_VectorLength = sqrt(r_VectorLength);


    size_t iterationsCnt = 0;
    const double EPSILON = 0.00001;

    // проверка на расхождение
    double prevEpsilonCheck = 0;
    size_t epsilonGrowCounter = 0;
    double epsilonCheck = r_VectorLength / b_VecLength;

    while (epsilonCheck >= EPSILON) {
        mulMatrixAndVector(matrixPart, N / size, z, mulMatrixAndVectorResult);
        MPI_Allgather(mulMatrixAndVectorResult, vectorPartCapacity, MPI_DOUBLE, Az, vectorPartCapacity, MPI_DOUBLE, MPI_COMM_WORLD);

        double r_r_partDotProductRes = dotProduct(r_vecPart, r_vecPart, vectorPartCapacity);
        double Az_Z_partDotProductRes = dotProduct(Az_vecPart, z_vecPart, vectorPartCapacity);
        //alpha = dotProduct(r, r) / dotProduct(Az, z);
        double alpha_part = r_r_partDotProductRes / Az_Z_partDotProductRes;
        mulVectorScalar(z, )
        //mulVectorScalar(z, alpha, alphaz);
        sumVector(x, alphaz, x);
        mulVectorScalar(Az, alpha, Az);
        memcpy(prev_r, r, sizeof(double) * N);
        subVector(r, Az, r);
        currDotProductRPart = dotProduct(r, r);
        beta = currDotProductRPart / prevDotProductRPart;
        prevDotProductRPart = currDotProductRPart;
        mulVectorScalar(z_vecPart, beta, betaz_vecPart, vectorPartCapacity);
        sumVector(r_vecPart, betaz_vecPart, z_vecPart, vectorPartCapacity);

        if(rank == 0){
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
        //MPI_Bcast(r, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Bcast(z, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //r_VectorLength = getVectorLength(r);
        r_VectorPartLength = getVectorPartSqr(r_vecPart, vectorPartCapacity);
        //double r_VectorLength = 0;
        MPI_Allreduce(&r_VectorPartLength, &r_VectorLength, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        r_VectorLength = sqrt(r_VectorLength);
        epsilonCheck = r_VectorLength / b_VecLength;
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
    //clock_t start = 0, end = 0;
    int size = 0, rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double* A = nullptr;
    double* b = nullptr;
    double* x = nullptr;

    if(rank == 0){
        A = new double[N * N];
        b = new double[N];
        x = new double[N];
    }

    MpiV2NonlinearConjugateGradient(A, b, x, rank, size);

    MPI_Finalize();

    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}
