#include <memory.h>
#include <cmath>
#include <random>
#include <iomanip>

#ifdef __unix__
#include <mpi.h>
#elif defined(_WIN32) || defined(WIN32)
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#endif

// divided by 1, 2, 4, 8, 16, 24
#define N (3840)
//#define N (8)
#define VAL_RANGE (50)

void initRandVector(double* vector) {
	for (size_t i = 0; i < N; ++i) {
		vector[i] = rand() % VAL_RANGE;
	}
}

void initVector(double* vector){
    for (size_t i = 0; i < N; ++i) {
        vector[i] = i+1;
    }
}
void initMatrix(double* matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i; j < N; ++j) {
            matrix[i * N + j] = i+1;
            matrix[j * N + i] = matrix[i * N + j];
        }
    }

    const int mainDiagonalWeighting = 1; // less number - run longer
    for (size_t i = 0; i < N; ++i) {
        matrix[i * N + i] += mainDiagonalWeighting;
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

// non-mpi version
void mulMatrixAndVector(const double* matrix, const double* vector, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = 0;
		for (size_t j = 0; j < N; ++j) {
			result[i] += vector[j] * matrix[i * N + j];
		}
	}
}

void mulMatrixAndVector(const double* matrix_part, int count, const double* vector, double* result) {
    memset(result, 0, count * sizeof(double));

    for (int i = 0; i < count; i++){
        for (int j = 0; j < N; j++) {
            result[i] += matrix_part[i * N + j] * vector[j];
        }
    }
}

void subVector(const double* first, const double* second, double* result) {
    for (size_t i = 0; i < N; ++i) {
        result[i] = first[i] - second[i];
    }
}

void subVector(const double* first, const double* second, double* result, const int& len) {
    for (size_t i = 0; i < len; ++i) {
        result[i] = first[i] - second[i];
    }
}

void sumVector(const double* first, const double* second, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = first[i] + second[i];
	}
}

void sumVector(const double* first, const double* second, double* result, const int& len) {
    for (size_t i = 0; i < len; ++i) {
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

double getVectorPartSqr(const double* vector, const int& len){
    double result = 0;
    for (size_t i = 0; i < len; ++i) {
        result += vector[i] * vector[i];
    }
    //return result/2;
    return result;
}

double getVectorLength(const double* vector, const int& len) {
    double result = 0;
    for (size_t i = 0; i < len; ++i) {
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

double dotProduct(const double* a, const double* b, const int& len) {
    double result = 0;
    for (size_t i = 0; i < len; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void mulVectorScalar(const double* vec, const double& scalar, double* result) {
	for (size_t i = 0; i < N; ++i) {
		result[i] = vec[i] * scalar;
	}
}

void mulVectorScalar(const double* vec, const double& scalar, double* result, const int& len) {
    for (size_t i = 0; i < len; ++i) {
        result[i] = vec[i] * scalar;
    }
}

void mulMatrixAndVectorParts(const int& size, const double* matrixPart, const double* vectorPart, double* result){
    double tmp[N] = {0};
    const int vecPartSize = N/size;
    const int matrixPartSize = N;

    for (int i = 0; i < vecPartSize; ++i) {
        for (int j = 0; j < matrixPartSize; ++j) {
            tmp[j] += matrixPart[i*N+j] * vectorPart[i];
        }
    }
    int* recvCounts = new int[size];
    for(int i = 0; i < size; ++i) recvCounts[i] = vecPartSize;
    MPI_Reduce_scatter(tmp, result, recvCounts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete[] recvCounts;
}

void printMatrix(const double* matrix, const size_t& len){
    for(size_t i = 0; i < len; ++i){
        for(size_t j = 0; j < len; ++j){
            std::cout << std::setw(2) << matrix[i * len + j] << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(const double* vec, const size_t& len){
    for(size_t i = 0; i < len; ++i){
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void initVector(double* vec, int val){
    for(size_t i= 0; i < N; ++i){
        vec[i] = val;
    }

}

void initMatrix(double* matrix, int val){
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            matrix[i*N + j] = val + (i == j);
        }
    }
}
