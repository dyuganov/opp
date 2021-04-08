#include <random>
#include <iostream>
#include <time.h>

#include "Matrix.h"

void MpiV2NonlinearConjugateGradient(double* A, double* b, double* x, int rank, int size){
    double startTime;
    if (rank == 0) startTime = MPI_Wtime();

    if (rank == 0) {
        /*initRandMatrix(A);
        initRandVector(b);
        initRandVector(x);*/
        initVector(b);
        initVector(x);
        initMatrix(A);
        /*initVector(b, 1);
        initVector(x, 1);
        initMatrix(A, 1);*/
    }

    if(rank == 0) std::cout << "hello0" << std::endl;

    const int matrixPartSize = N * N / size;
    const int vectorPartSize = N / size;
    double alpha = 0, beta = 0;
    double* A_matrixPart = new double[matrixPartSize];
    double* Az_vecPart = new double[vectorPartSize];
    double* alphaz_vecPart = new double[vectorPartSize];
    double* betaz_vecPart = new double[vectorPartSize];
    double* b_vecPart = new double[vectorPartSize];
    double* x_vecPart = new double[vectorPartSize];
    double* z_vecPart = new double[vectorPartSize];
    double* Ax_vecPart = new double[vectorPartSize];
    double* Az_tmpVecPart = new double[vectorPartSize];
    double* Ax_tmpVecPart = new double[vectorPartSize];
    double* r_vecPart = new double[vectorPartSize];
    double* r_prevVecPart = new double[vectorPartSize];

    double r_prevDotProduct = 0;
    double r_prevDotProductPart = 0;
    double r_dotProduct = 0;
    double r_dotProductPart = 0;
    double r_vecLength = 0;
    double r_vecPartLength = 0;
    double b_vecLength = 0;
    if(rank == 0){
        b_vecLength = getVectorLength(b);
    }
    MPI_Bcast(&b_vecLength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // один раз посчитали длину b раздали всем

    if(rank == 0) std::cout << "hello1" << std::endl;
    for(int i = 0; i < vectorPartSize; ++i){
        for(int j = 0; j < vectorPartSize; ++j){
            A_matrixPart[i * vectorPartSize + j] = 0;
        }
    }
    for(int i = 0; i < vectorPartSize; ++i){
        Ax_vecPart[i] = 0;
        z_vecPart[i] = 0;
        betaz_vecPart[i] = 0;
        alphaz_vecPart[i] = 0;
        Az_vecPart[i] = 0;
        b_vecPart[i] = 0;
        x_vecPart[i] = 0;
        Ax_tmpVecPart[i] = 0;
        Az_tmpVecPart[i] = 0;
        r_vecPart[i] = 0;
        r_prevVecPart[i] = 0;
    }
    if(rank == 0) std::cout << "hello2" << std::endl;
    MPI_Scatter(x, vectorPartSize, MPI_DOUBLE, x_vecPart, vectorPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, vectorPartSize, MPI_DOUBLE, b_vecPart, vectorPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(A, matrixPartSize, MPI_DOUBLE, A_matrixPart, matrixPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(rank == 0) std::cout << "hello3" << std::endl;
    //mulMatrixAndVectorParts(A_matrixPart, x_vecPart, Ax_vecPart, size, rank); // Ax^0
    //subVector(b_vecPart, Ax_vecPart, r_vecPart, vectorPartSize); // r^0 = b - Ax^0

    //mulMatrixAndVectorParts(A_matrixPart, x_vecPart, Ax_vecPart, size, rank); // Ax^0
    //memset(Ax_vecPart, 0, vectorPartSize);
    for(size_t i = 0; i < vectorPartSize; ++i) Ax_vecPart[i] = 0;
    int vecPartsSize[vectorPartSize];
    if(rank == 0) std::cout << "hello4" << std::endl;
    for(size_t i = 0; i < vectorPartSize; ++i) vecPartsSize[i] = vectorPartSize;
    double* tmp = new double[vectorPartSize];
    if(rank == 0) std::cout << "hello5" << std::endl;
    mulMatrixPartAndVectorPart(A_matrixPart, x_vecPart, Ax_vecPart, size, rank);
    if(rank == 0) std::cout << "hello6" << std::endl;
    if(rank == 0) std::cout << "vecPartsSize "  << vectorPartSize << std::endl;
    //MPI_Reduce_scatter(Ax_vecPart, tmp, vecPartsSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // problem here
    double* Ax = new double[N];
    MPI_Reduce(Ax_vecPart, Ax, vectorPartSize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Scatter(Ax, vectorPartSize, MPI_DOUBLE, Ax_vecPart, vectorPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(rank == 0) std::cout << "hello7" << std::endl;
    //memcpy(tmp, Ax_vecPart, vectorPartSize);
    for(size_t i = 0; i < vectorPartSize; ++i) Ax_vecPart[i] = tmp[i];

    /*double* matrix = new double[N*N];
    MPI_Allgather(Ax_vecPart, N*N/size, MPI_DOUBLE, matrix, N*N, MPI_DOUBLE, MPI_COMM_WORLD);
    if (rank == 0) std::cout << "matrix: ";
    if (rank == 0) for(int i = 0; i < matrixPartSize; ++i) std::cout << Ax_vecPart[i] << ' ';
    if (rank == 0) std::cout << std::endl;*/
    std::cout << "hello";

    subVector(b_vecPart, Ax_vecPart, r_vecPart, vectorPartSize); // r^0 = b - Ax^0

    for(int i = 0; i < vectorPartSize; ++i) {
        z_vecPart[i] = r_vecPart[i]; // z^0 = r^0
        r_prevVecPart[i] = r_vecPart[i];
    }

    // (r^n, r^n)
    r_prevDotProductPart = dotProduct(r_vecPart, r_vecPart, vectorPartSize);
    MPI_Allreduce(&r_prevDotProductPart, &r_prevDotProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //if(rank == 0) std::cout << "Dot 1: " << r_prevDotProduct << std::endl;

    // длина r для первой проверки
    r_vecPartLength = getVectorPartSqr(r_vecPart, vectorPartSize);
    MPI_Allreduce(&r_vecPartLength, &r_vecLength, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    r_vecLength = sqrt(r_vecLength);

    //if(rank == 0) std::cout << "r_vecLength: " << r_vecLength << std::endl;
    //if(rank == 0) std::cout << "b_vecLength: " << b_vecLength << std::endl;

    size_t iterationsCnt = 0;
    const double EPSILON = 0.00001;

    // проверка на расхождение
    double prevEpsilonCheck = 0;
    size_t epsilonGrowCounter = 0;
    double epsilonCheck = r_vecLength / b_vecLength; // first EPSILON check

    double Az_z_dotProductRes = 0, Az_z_partDotProductRes = 0;

    while(epsilonCheck >= EPSILON){
        //Az^n
        //mulMatrixAndVectorParts(A_matrixPart, z_vecPart, Az_vecPart, size, rank);
        for(size_t i = 0; i < vectorPartSize; ++i) Az_vecPart[i] = 0;

        mulMatrixPartAndVectorPart(A_matrixPart, z_vecPart, Az_vecPart, size, rank);
        //MPI_Reduce_scatter(Az_vecPart, tmp, vecPartsSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double* Az = new double[N];
        MPI_Reduce(Az_vecPart, Az, vectorPartSize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Scatter(Az, vectorPartSize, MPI_DOUBLE, Az_vecPart, vectorPartSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //memcpy(tmp, Az_vecPart, vectorPartSize);
        //for(size_t i = 0; i < vectorPartSize; ++i) Az_vecPart[i] = tmp[i];

        // (Az^n, z^n)
        Az_z_dotProductRes = 0;
        Az_z_partDotProductRes = dotProduct(Az_vecPart, z_vecPart, vectorPartSize);
        MPI_Allreduce(&Az_z_partDotProductRes, &Az_z_dotProductRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //if(rank == 0) std::cout << "Az_z_dotProductRes " << Az_z_dotProductRes << std::endl;

        // alpha^n+1 = (r^n, r^n) / (Az^n, z^n)
        //alpha = r_prevDotProduct / Az_z_dotProductRes;
        alpha = r_dotProduct / Az_z_dotProductRes;
        //if(rank == 0) std::cout << "r_prevDotProduct " << r_prevDotProduct << std::endl;
        //if(rank == 0) std::cout << "Az_z_dotProductRes " << Az_z_dotProductRes << std::endl;
        //if(rank == 0) std::cout << "alpha " << alpha << std::endl;

        //x^n+1 = x^n + alpha^n+1 * z^n
        sumVector(x_vecPart, alphaz_vecPart, x_vecPart, vectorPartSize);

        //alpha^n+1 * z^n
        mulVectorScalar(Az_vecPart, alpha, Az_vecPart, vectorPartSize);
        //mulVectorScalar(z_vecPart, alpha, alphaz_vecPart,vectorPartSize);

        // r^n+1 = r^n - alpha^n+1 * z^n
        subVector(r_vecPart, Az_vecPart, r_vecPart, vectorPartSize);
        //subVector(r_vecPart, alphaz_vecPart, r_vecPart, vectorPartSize);

        // (r^n+1, r^n+1)
        r_dotProductPart = dotProduct(r_vecPart, r_vecPart, vectorPartSize);
        MPI_Allreduce(&r_dotProductPart, &r_dotProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //if(rank == 0) std::cout << "r_dotProduct " << r_dotProduct << std::endl;

        // beta^n+1 = (r^n+1, r^n+1) / (r^n, r^n)
        beta = r_dotProduct / r_prevDotProduct;
        r_prevDotProduct = r_dotProduct;
        //if(rank == 0) std::cout << "beta " << beta << std::endl;

        // beta^n+1 * z^n
        mulVectorScalar(z_vecPart, beta, betaz_vecPart, vectorPartSize);

        // z^n+1 = r^n+1 + beta^n+1 * z^n
        sumVector(r_vecPart, betaz_vecPart, z_vecPart, vectorPartSize);

        // ||r^n||
        r_vecPartLength = getVectorPartSqr(r_vecPart, vectorPartSize);
        MPI_Allreduce(&r_vecPartLength, &r_vecLength, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        r_vecLength = sqrt(r_vecLength);
        //if(rank == 0) std::cout << "r_vecLength " << r_vecLength << std::endl;
        //if(rank == 0) std::cout << "b_vecLength " << r_vecLength << std::endl;
        epsilonCheck = r_vecLength / b_vecLength;

        // resolve check + iterations counter
        if(rank == 0){
            std::cout << "epsilonCheck " << epsilonCheck << std::endl;
            if(epsilonCheck >= prevEpsilonCheck){
                ++epsilonGrowCounter;
                prevEpsilonCheck = epsilonCheck;
            }
            if(epsilonGrowCounter > 5){
                perror("Can't resolve the matrix.");
                std::cout << "Can't resolve the matrix." << std::endl;
                std::cout << "Iterations: " << iterationsCnt << std::endl;
                std::cout << "TIME: " << MPI_Wtime() - startTime << std::endl;
                exit(1);
            }
            ++iterationsCnt;
        }
    }

    if(rank == 0) {
        std::cout << "Iterations: " << iterationsCnt << std::endl;
        std::cout << "TIME: " << MPI_Wtime() - startTime << std::endl;
    }

    delete[] A_matrixPart;
    delete[] Az_vecPart;
    delete[] alphaz_vecPart;
    delete[] betaz_vecPart;
    delete[] b_vecPart;
    delete[] x_vecPart;
    delete[] z_vecPart;
    delete[] Ax_vecPart;
    delete[] Az_tmpVecPart;
    delete[] Ax_tmpVecPart;
    delete[] r_vecPart;
    delete[] r_prevVecPart;
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


/*    const int vecLen = 8;
    double* vec1 = new double[vecLen];
    for (int i = 0; i < 8; ++i) vec1[i] = i + 1;
    double* vec2 = new double[vecLen];
    for (int i = 0; i < 8; ++i) vec2[i] = 0;

    double len = 0;
    double vecPartLen = getVectorPartSqr(vec1, vecLen);
    MPI_Allreduce(&vecPartLen, &len, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std::cout << "vec: ";
    for (int i = 0; i < 8; ++i) std::cout << vec1[i] << ' ';
    len = sqrt(len)/2;
    std::cout << std::endl<< "len: " << len/2 << std::endl;*/


    MPI_Finalize();

    if(rank == 0){
        delete[] A;
        delete[] b;
        delete[] x;
    }

    return 0;
}