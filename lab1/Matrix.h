#pragma once

// divided by 1, 2, 4, 8, 16, 24
#define N (3840)

void initRandVector(double* vector) {
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

	const int mainDiagonalWeighting = 400; // less number - run longer
	for (size_t i = 0; i < N; ++i) {
		matrix[i * N + i] += mainDiagonalWeighting;
	}
}

// non-mpi version
void mulMatrixAndVector(const double* matrix, const double* vector, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = 0;
		for (size_t j = 0; j < N; ++j) {
			result[i] += vector[j] * matrix[i * N + j];
		}
	}
}
// mpi version
void mulMatrixAndVector(const double* matrix, const double* vector, double* result, const int& rank) {
    for (size_t i = 0; i < N; ++i) {
        result[i] = 0;
        for (size_t j = 0; j < N; ++j) {
            result[i] += vector[j] * matrix[i * N + j];
        }
    }
}

void sumVector(const double* first, const double* second, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = first[i] + second[i];
	}
}

void subVector(const double* first, const double* second, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = first[i] - second[i];
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

void mulVectorScalar(const double* matrix, const double& scalar, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = matrix[i] * scalar;
	}
}