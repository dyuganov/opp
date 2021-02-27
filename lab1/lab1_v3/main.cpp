#include<mpi.h> // Подключение библиотеки MPI
#include<stdio.h>
#include <random>

#define N (3000)
#define EPSILON (0.00001)

void initRandVector(double* vector){
	for (size_t i = 0; i < N; ++i) {
		vector[i] = rand() % 10;
	}
}

void initRandMatrix(double* matrix) {
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = i; j < N; ++j) {
			matrix[i * N + j] = rand() % 10;
			matrix[j * N + i] = matrix[i * N + j];
		}
	}

	const int mainDiagonalWeighting = 100;
	for (size_t i = 0; i < N; ++i) {
		matrix[i * N + i] += mainDiagonalWeighting;
	}
}

void mulMatrixAndVector(double* matrix, double* vector, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = 0;
		for (size_t j = 0; j < N; ++j) {
			result[i] += vector[j] * matrix[i * N + j];
		}
	}
}

void sumVector(double* first, double* second, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = first[i] + second[i];
	}
}

void subVector(double* first, double* second, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = first[i] - second[i];
	}
}

double getVectorLength(double* vector) {

}

//void subMatrixAndVector(double* matrix, double* vector, double* result) {
//
//}
//
//void sumMatrixAndVector(double* matrix, double* vector, double* result) {
//
//}

double doProduct(double* a, double* b) {
	double dp = 0;
	for (size_t i = 0; i < N; ++i) {
		dp += a[i] * b[i];
	}
	return dp;
}

void mulScalar(double* matrix, const double& scalar, double* result) {

}

double* nonlinearConjugateGradient(double* A, double* b) {
	double* r = new double[N];
	double* z = new double[N];
	double* result = new double[N];

	initRandMatrix(A);
	initRandVector(b);
	initRandVector(result);

	double* tmp = new double[N];
	mulMatrixAndVector(A, result, tmp);
	subVector(b, tmp, r);
	delete[] tmp;

	memcpy(z, r, sizeof(double) * N);


	double alpha = 0;
	double beta = 0;
	double* Az = new double[N];
	double* alphaz = new double[N];
	double* betaz = new double[N];
	double* prev_r = new double[N];
	
	const double vectorBLength = getVectorLength(b);
	while (getVectorLength(r) / vectorBLength >= EPSILON) {

		mulMatrixAndVector(A, z, Az);
		alpha = doProduct(r, r) / doProduct(Az, z);

		mulScalar(z, alpha, alphaz);

		sumVector(result, alphaz, result);
		mulScalar(Az, alpha, Az);

		memcpy(prev_r, r, sizeof(double) * N);
		subVector(r, Az, r);
		beta = doProduct(r, r) / doProduct(prev_r, prev_r);

		mulScalar(z, beta, betaz);
		sumVector(r, betaz, z);
	}

	delete[] alphaz;
	delete[] Az;
	delete[] betaz;
	delete[] prev_r;
	delete[] r;
	delete[] z;

	return result;
}

void noMPI() {
	double* A = new double[N*N];
	double* b = new double[N];
	double* x = nullptr;

	x = nonlinearConjugateGradient(A, b);

	delete[] A;
	delete[] b;
	delete[] x;
}


int main(int argc, char* argv[]) {

	int size, rank;
	MPI_Init(&argc, &argv); // Инициализация MPI

	double* A = nullptr;
	double* b = nullptr;
	double* x = nullptr;


	MPI_Finalize(); // Завершение работы MPI

	return 0;
}

