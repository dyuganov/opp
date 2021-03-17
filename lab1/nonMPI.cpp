#include <random>
#include <iostream>
#include <time.h>
#include <memory.h>
#include <cmath>

// divided by 1, 2, 4, 8, 16, 24
#define N (3840)
#define VAL_RANGE (100)

void initRandVector(double* vector) {
    for (size_t i = 0; i < N; ++i) {
        vector[i] = rand() % VAL_RANGE;
    }
}

void initRandMatrix(double* matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i; j < N; ++j) {
            matrix[i * N + j] = rand() % VAL_RANGE;
            matrix[j * N + i] = matrix[i * N + j];
        }
    }

    const int mainDiagonalWeighting = 350; // less number - run longer
    for (size_t i = 0; i < N; ++i) {
        matrix[i * N + i] += mainDiagonalWeighting;
    }
}

void mulMatrixAndVector(const double* matrix, const double* vector, double* result) {
    for (size_t i = 0; i < N; ++i) {
        result[i] = 0;
        for (size_t j = 0; j < N; ++j) {
            result[i] += vector[j] * matrix[i * N + j];
        }
    }
}

void subVector(const double* first, const double* second, double* result) {
    for (size_t i = 0; i < N; ++i) {
        result[i] = first[i] - second[i];
    }
}

void sumVector(const double* first, const double* second, double* result) {
    for (size_t i = 0; i < N; ++i) {
        result[i] = first[i] + second[i];
    }
}

double getVectorLength(const double* vector) {
    double result = 0;
    for (size_t i = 0; i < N; ++i) {
        result += vector[i] * vector[i];
    }
    return sqrt(result);
}

double dotProduct(const double* a, const double* b) {
    double result = 0;
    for (size_t i = 0; i < N; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void mulVectorScalar(const double* vec, const double& scalar, double* result) {
    for (size_t i = 0; i < N; ++i) {
        result[i] = vec[i] * scalar;
    }
}

double* nonMpiNonlinearConjugateGradient(double* A, double* b, double* x) {
	double* r = new double[N];
	double* z = new double[N];
	double alpha = 0, beta = 0;
	double* Az = new double[N];
	double* alphaz = new double[N];
	double* betaz = new double[N];
	double* prev_r = new double[N];

	initRandMatrix(A);
	initRandVector(b);
	initRandVector(x);

	double* tmp = new double[N];
	mulMatrixAndVector(A, x, tmp);
	subVector(b, tmp, r);
	delete[] tmp;

	memcpy(z, r, sizeof(double) * N);

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

	double* A = new double[N * N];
	double* b = new double[N];
	double* x = new double[N];

	clock_t start, end;

	start = clock();
	nonMpiNonlinearConjugateGradient(A, b, x);
	end = clock();

	std::cout << "Time: " << (double)((end - start) / CLOCKS_PER_SEC) << " sec." << std::endl;

	delete[] A;
	delete[] b;
	delete[] x;

	return 0;
}
