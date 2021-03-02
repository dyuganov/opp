#include <random>
#include <iostream>
#include <time.h>

//#include <mpi.h> // for cluster
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h" // for local use

#include "Matrix.h"

// vectors duplicate in each process
double* MpiV1NonlinearConjugateGradient(double* A, double* b, double* x, int rank, int size) {
	if (rank == 0) {
		initRandMatrix(A);
		initRandVector(b);
		initRandVector(x);
	}

	
	int matrixPartCapacity = N * N / size;
	//int vectorPartCapacity = N / size;
	auto* matrixPart = new double[matrixPartCapacity];
	
	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD); // раздали переменные
	MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// раздали матрицу по кускам
	MPI_Scatter(A, matrixPartCapacity, MPI_DOUBLE, matrixPart, matrixPartCapacity, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	auto* r = new double[N];
	auto* z = new double[N];

	auto* tmp = new double[N];
	mulMatrixAndVector(A, x, tmp);
	subVector(b, tmp, r);
	delete[] tmp;

	memcpy(z, r, sizeof(double) * N);

	double alpha = 0, beta = 0;
	auto* Az = new double[N];
	auto* alphaz = new double[N];
	auto* betaz = new double[N];
	auto* prev_r = new double[N];

	double prevDoProductR = dotProduct(r, r);
	double currDoProductR;
	const double bVectorLength = getVectorLength(b);
	double rVectorLength = getVectorLength(r);
	size_t cnt = 0;
	const double epsilon = 0.00001;
	while (rVectorLength / bVectorLength >= epsilon) {
		mulMatrixAndVector(A, z, Az);
		alpha = dotProduct(r, r) / dotProduct(Az, z);

		mulVectorScalar(z, alpha, alphaz);

		sumVector(x, alphaz, x);
		mulVectorScalar(Az, alpha, Az);

		memcpy(prev_r, r, sizeof(double) * N);

		subVector(r, Az, r);

		currDoProductR = dotProduct(r, r);
		beta = currDoProductR / prevDoProductR;
		prevDoProductR = currDoProductR;

		mulVectorScalar(z, beta, betaz);
		sumVector(r, betaz, z);

		rVectorLength = getVectorLength(r);
		++cnt;
	}

	std::cout << "Iterations: " << cnt << std::endl;

	delete[] alphaz;
	delete[] Az;
	delete[] betaz;
	delete[] prev_r;
	delete[] r;
	delete[] z;

	return x;
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


	if(rank == 0) start = clock();

	MpiV1NonlinearConjugateGradient(A, b, x, rank, size);

	if (rank == 0) end = clock();
	if (rank == 0) std::cout << "Time: " << ((double)(end - start) / CLOCKS_PER_SEC) << " sec." << std::endl;

	delete[] A;
	delete[] b;
	delete[] x;

	MPI_Finalize();

	return 0;
}
