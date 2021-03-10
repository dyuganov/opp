#include <random>
#include <iostream>
#include <time.h>

#include "Matrix.h"

// vectors duplicate in each process
void MpiV2NonlinearConjugateGradient(double* A, double* b, double* x, int rank, int size) {
    double startTime;
    if (rank == 0) startTime = MPI_Wtime();

    if (rank == 0) {
        initRandMatrix(A);
        initRandVector(b);
        initRandVector(x);
    }

    double alpha = 0, beta = 0;

    const int matrixPartSize = N * N / size;
    const int vectorPartSize = N / size;
    auto* A_matrixPart = new double[matrixPartSize];

    auto* Az_vecPart = new double[vectorPartSize];
    auto* alphaz_vecPart = new double[vectorPartSize];
    auto* betaz_vecPart = new double[vectorPartSize];
    auto* b_vecPart = new double[vectorPartSize];
    auto* x_vecPart = new double[vectorPartSize];
    auto* z_vecPart = new double[vectorPartSize];
    auto* Ax_vecPart = new double[vectorPartSize];

    MPI_Scatter(x, vectorPartSize, MPI_DOUBLE, x_vecPart, vectorPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, vectorPartSize, MPI_DOUBLE, b_vecPart, vectorPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(A, matrixPartSize, MPI_DOUBLE, A_matrixPart, matrixPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    mulMatrixAndVectorParts(A_matrixPart, x_vecPart, Ax_vecPart, size, rank);

    auto* r_vecPart = new double[vectorPartSize];
    auto* r_prevVecPart = new double[vectorPartSize];
    subVector(b_vecPart, Ax_vecPart, r_vecPart, vectorPartSize);

    //MPI_Scatter(r, vectorPartSize, MPI_DOUBLE, r_vecPart, vectorPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Scatter(z, vectorPartSize, MPI_DOUBLE, z_vecPart, vectorPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    memcpy(z_vecPart, r_vecPart, sizeof(double) * vectorPartSize);
    // подготовили начало, нуллевые значения

    double b_vecLength = 0;
    if(rank == 0){
        b_vecLength = getVectorLength(b);
    }
    MPI_Bcast(&b_vecLength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // один раз посчитали длину b раздали всем

    // (r, r)
    double r_prevDotProductPart = dotProduct(r_vecPart, r_vecPart, vectorPartSize);
    double r_prevDotProduct = 0;
    MPI_Allreduce(&r_prevDotProductPart, &r_prevDotProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double r_currDotProductPart = 0, r_currDotProduct = 0;

    // длина r для первой проверки
    double r_vecLength = 0;
    double r_vecPartLength = getVectorPartSqr(r_vecPart, vectorPartSize);
    MPI_Allreduce(&r_vecPartLength, &r_vecLength, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    r_vecLength = sqrt(r_vecLength);

    size_t iterationsCnt = 0;
    const double EPSILON = 0.00001;

    // проверка на расхождение
    double prevEpsilonCheck = 0;
    size_t epsilonGrowCounter = 0;
    double epsilonCheck = r_vecLength / b_vecLength; // первая проверка на EPSILON

    while (epsilonCheck >= EPSILON) {
        mulMatrixAndVectorParts(A_matrixPart, z_vecPart, Az_vecPart, size, rank);
        // получили куски вектора Az в каждом потоке

        //r_prevDotProductPart = dotProduct(r_vecPart, r_vecPart, vectorPartSize);
        //MPI_Allreduce(&r_prevDotProduct, &r_prevDotProductPart, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double Az_z_dotProductRes = 0;
        double Az_z_partDotProductRes = dotProduct(Az_vecPart, z_vecPart, vectorPartSize);
        MPI_Allreduce(&Az_z_partDotProductRes, &Az_z_dotProductRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //alpha = r_prevDotProductPart / Az_z_dotProductRes; // собрали альфу по частям, нашли
        alpha = r_currDotProduct / Az_z_dotProductRes; // собрали альфу по частям, нашли


        mulVectorScalar(z_vecPart, alpha, alphaz_vecPart, vectorPartSize); // альфа и z перемножили, получили alphaz по кускам в потоках
        sumVector(x_vecPart, alphaz_vecPart, x_vecPart, vectorPartSize);
        mulVectorScalar(Az_vecPart, alpha, Az_vecPart, vectorPartSize);
        memcpy(r_prevVecPart, r_vecPart, sizeof(double) * vectorPartSize);
        subVector(r_prevVecPart, Az_vecPart, r_vecPart, vectorPartSize);

        r_currDotProductPart = dotProduct(r_vecPart, r_vecPart, vectorPartSize);
        MPI_Allreduce(&r_currDotProductPart, &r_currDotProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta = r_currDotProduct / r_prevDotProduct;
        r_prevDotProduct = r_currDotProduct;
        r_prevDotProductPart = r_currDotProduct;

        mulVectorScalar(z_vecPart, beta, betaz_vecPart, vectorPartSize);

        sumVector(r_vecPart, betaz_vecPart, z_vecPart, vectorPartSize);

        if(rank == 0){
            if(epsilonCheck > prevEpsilonCheck){
                ++epsilonGrowCounter;
                prevEpsilonCheck = epsilonCheck;
            }
            //std::cout << "r_length: " << r_vecLength << std::endl;
            if(epsilonGrowCounter > 5){
                perror("Can't resolve the matrix.");
                std::cout << "Can't resolve the matrix." << std::endl;
                std::cout << "Iterations: " << iterationsCnt << std::endl;
                std::cout << "TIME: " << MPI_Wtime() - startTime << std::endl;
                exit(1);
            }
            ++iterationsCnt;
        }

        r_vecPartLength = getVectorPartSqr(r_vecPart, vectorPartSize);
        MPI_Allreduce(&r_vecPartLength, &r_vecLength, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        r_vecLength = sqrt(r_vecLength);

        epsilonCheck = r_vecLength / b_vecLength;
    }

    if(rank == 0) {
        std::cout << "Iterations: " << iterationsCnt << std::endl;
        std::cout << "TIME: " << MPI_Wtime() - startTime << std::endl;
    }

    delete[] r_vecPart;
    delete[] r_prevVecPart;
    delete[] A_matrixPart;
    delete[] Az_vecPart;
    delete[] alphaz_vecPart;
    delete[] betaz_vecPart;
    delete[] b_vecPart;
    delete[] x_vecPart;
    delete[] z_vecPart;
    delete[] Ax_vecPart;
}

int main(int argc, char* argv[]) {
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
