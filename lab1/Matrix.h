#pragma once

#include <memory.h>
#include <cmath>
#include <random>

//#include <mpi.h> // for cluster
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h" // for local use

// divided by 1, 2, 4, 8, 16, 24
#define N (3840)
#define VAL_RANGE (20)

void initRandVector(double* vector) {
	for (size_t i = 0; i < N; ++i) {
		vector[i] = rand() % VAL_RANGE;
	}
}

void initVector(double* vector, const double& val){
    for (size_t i = 0; i < N; ++i) {
        vector[i] = val;
    }
}

void initEMatrix(double* matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i; j < N; ++j) {
            matrix[i * N + j] = 1;
            matrix[j * N + i] = matrix[i * N + j];
        }
    }

    for (size_t i = 0; i < N; ++i) {
        matrix[i * N + i] = 2;
    }
}

void initRandMatrix(double* matrix) {
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = i; j < N; ++j) {
			matrix[i * N + j] = rand() % VAL_RANGE;
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
/*void mulMatrixAndVector(const double* matrix, const double* vector, double* result, const int& rank) {
    for (size_t i = 0; i < N; ++i) {
        result[i] = 0;
        for (size_t j = 0; j < N; ++j) {
            result[i] += vector[j] * matrix[i * N + j];
        }
    }
}*/

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

/*void mulMatrixAndVectorParts(const double* matrixPart, double* vectorPart, double* result, const int& size, const int& rank) {
    int rows = N / size;
    const int begin = rows * rank;
    const int length = rows;
    const int send_len = (int) N / size;

    for (int shift = 0; shift < size; ++shift) {
        for (int i = 0; i < rows; ++i) {
            for (int j = begin; j < begin + length; ++j) {
                result[i] += matrixPart[i * N + j] * vectorPart[j - begin]; // BEGIN not changes!!!
            }
        }
        int send_id = (rank + 1) % size;
        int recv_id = (rank + size - 1) % size;

        //send curr x to another proc
        const int tag = 42;
        MPI_Sendrecv_replace(vectorPart, send_len, MPI_DOUBLE, send_id, tag, recv_id, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}*/

void mulMatrixAndVectorParts(const double* matrixPart, double* vectorPart, double* result, const int& size, const int& rank) {
    int rows = N / size;
    int length = (int)rows;
    int begin = length * rank;
    const int send_len = (int) N / size;
    //int sum = 0;

    for (int shift = 0; shift < size; ++shift) {
        for (int i = 0; i < rows; ++i) {
            for (int j = begin; j < begin + length; ++j) {
                result[i] += matrixPart[i * N + j] * vectorPart[j - begin];
            }
        }

        int send_id = (rank + 1) % size;
        int recv_id = (rank + size - 1) % size;

        int pos = (rank + size - shift - 1) % size;
        begin = length * pos;
        length = length % N;

        const int tag = 42;
        MPI_Sendrecv_replace(vectorPart, send_len, MPI_DOUBLE, send_id, tag, recv_id, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/*void mulMatrixAndVectorParts(const double* matrixPart, double* vectorPart, double* result, const int& size, const int& rank) {
    int rows = N / size;
    int length =  N / size;
    int begin = rows * rank;
    const int send_len = (int) N / size;

    for (int shift = 0; shift < size; ++shift) {
        for (int i = 0; i < rows; ++i) {
            for (int j = begin; j < begin + length; ++j) {
                result[i] += matrixPart[i * N + j] * vectorPart[j - begin];
            }
        }

        int pos = (rank + size - shift - 1) % size;
        begin = pos * length;
        length = length % N;

        *//*begin = pos * length;
        length = length % N;*//*

        int send_id = (rank + 1) % size;
        int recv_id = (rank + size - 1) % size;


        const int tag = 42;
        MPI_Sendrecv_replace(vectorPart, send_len, MPI_DOUBLE, send_id, tag, recv_id, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}*/

void mulColAndVecPart(double* matrixPart, double* vecPart, double* resultVecPart, int vecPartSize){
    for (int i = 0, strCnt = 0; i < vecPartSize; ++i){
        if (i % N == 0){
            ++strCnt;
        }
        resultVecPart[i % N] += matrixPart[i] * vecPart[strCnt];
    }
}

