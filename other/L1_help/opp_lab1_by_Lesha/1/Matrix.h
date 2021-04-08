#include "Constants.h"
#include <cstring>

int* getAllLengths(int size){ //функция составления массива длин (в строках) частей матрицы в зависимости от количества процессов
    int* result = (int*) malloc(sizeof(int) * size);
    int add = N % size;
    for (size_t i = 0; i < size; i++)
        result[size-i-1] = N / size + (add-- > 0 ? 1 : 0);
    return result;
}

int* prepareDispls(int size, int* lengths) {//функция составления массива отступов (в строках) до частей матрицы в зависимости от количества процессов
    int* result = (int*) malloc(sizeof(int) * size);
    result[0] = 0;
    int sum = 0;
    for (int i = 1; i < size; i++) {
        sum += lengths[i - 1];
        result[i] = sum;
    }
    return result;
}

void matrixInit(double* matrix) {//базовое заполнение матрицы по условию задачи
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            matrix[i*N + j] = i == j ? 2 : 1;
}

void otherInit(double* b, double* x){//базовое заполнение векторов b и x
    for(int i = 0; i < N; i++){
        b[i] = N + 1;
        x[i] = 0;
    }
}

double getVecLen(double* vector){// функция подсчета длины вектора
    double sum = 0;
    for (int i = 0; i < N; i++)
        sum += vector[i] * vector[i];
    return sqrt(sum);
}

void vecSub (double * vec1 , double* vec2, double* res){// функция разности векторов
    for (int i = 0; i < N; i++) {
        res[i] = vec1[i] - vec2[i];
    }
}

void vecMulScalar (double scalar, double* vec, double* res){// функция умножения вектора на скаляр
    for (int i = 0; i < N; i++){
        res[i] = scalar * vec[i];
    }
}

bool isEnd (double sum, double bLength){//функция проверки на конец итерирования программы
    return ((sum) / (bLength)) < e;
}

void printVector(double* vec) {//функция печати вектора
    for (int i = 1; i <= N; i++) {
        printf("%.5f ", vec[i - 1]);
        if (i % N == 0) printf("\n");
    }
}

void mulMatrixPart(double* matrix_part, int count, double* vector, double* result) { //функция умножения части матрицы на вектор
    memset(result, 0, count * sizeof(double));

    for (int i = 0; i < count; i++)
        for (int j = 0; j < N; j++)
            result[i] += matrix_part[i*N + j] * vector[j];
}

void freeAll(int* f, int* g, double* a, double* b, double*c,
             double* d, double* j, double* k,  int* h, int* i){//а что выглядит хайпово
    free(a);
    free(b);
    free(c);
    free(d);
    free(f);
    free(g);
    free(h);
    free(i);
    free(j);
    free(k);
}