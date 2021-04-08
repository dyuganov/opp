#include "Constants.h"

void getAllLengths(int size, int* result){ //функция составления массива длин (в строках) частей матрицы в зависимости от количества процессов
    int add = N % size;
    for (int i = 0; i < size; i++)
        result[size-i-1] = N / size + (add-- > 0 ? 1 : 0);
}

void prepareDispls(int size, int* lengths, int* result) { //функция составления массива отступов (в строках) до частей матрицы в зависимости от количества процессов
    result[0] = 0;
    int sum = 0;
    for (int i = 1; i < size; i++) {
        sum += lengths[i - 1];
        result[i] = sum;
    }
}

void matrixInit(double* matrix) {//базовое заполнение матрицы по условию задачи
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            matrix[i*N + j] = i == j ? 2 : 1;
}

void initBX(double* full_x, double* full_b){//базовое заполнение векторов b и x
    for(int i = 0; i < N; i++){
        full_b[i] = N + 1;
        full_x[i] = 0;
    }
}

double getVecLen(double* vector, int count){// функция подсчета длины вектора
    double sum = 0;
    for (int i = 0; i < count; i++)
        sum += vector[i] * vector[i];
    return sqrt(sum);
}

void vecSub (double * vec1 , double* vec2, double* res, int size){ // функция разности векторов
    for (int i = 0; i < size; i++) {
        res[i] = vec1[i] - vec2[i];
    }
}

void vecMulScalar (double* vec, double scalar, double* res, int size){ // функция умножения вектора на скаляр
    for (int i = 0; i < size; i++){
        res[i] = scalar * vec[i];
    }
}

bool isEnd (double sum, double bLength){ //функция проверки на конец итерирования программы
    return (sqrt(sum) / bLength) < e;
}

void printVector(double* vec, int count) {  //функция печати вектора
    for (int i = 1; i <= count; i++) {
        printf("%.5f ", vec[i - 1]);
        if (i % N == 0) printf("\n");
    }
}

void mulMatrixPart(double* matrix, double* vector, double* result, int* lengths, int* displs, int size, int rank) {
    int rows = lengths[rank]; //кол-во строк
    int begin = displs[rank]; //начало в зависимости от отступа
    int length = rows; //
    int send_len = (int) ceil((double) N / size); // a.k.a extended

    for (int shift = 0; shift < size; shift++) { // пересчет Ax для пересылки на каждый процесс, чтобы в конце сумма на
        for (int i = 0; i < rows; i++)           // каждом процессе была одинаковая, иначе в силу разницы в элементах матрицы разделять вектор бы не получилось
            for (int j = begin; j < begin + length; j++)
                result[i] += matrix[i*N + j] * vector[j - begin];

        int pos = (rank + size - shift - 1) % size;
        begin = displs[pos];
        length = lengths[pos] % N;
        int send_id = (rank + 1) % size;
        int recv_id = (rank + size - 1) % size;

        MPI_Sendrecv_replace(vector, send_len, MPI_DOUBLE, send_id,5,
                             recv_id, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//пересылаем текущий x а другой процесс
    }

}

void freeAll( double* a, double* b, double*c ,double* d, int* f,
              int* g, int* h, int* i, double* k, double* l, double* y){ //а что выглядит хайпово
    free(a);
    free(b);
    free(c);
    free(d);
    free(f);
    free(g);
    free(h);
    free(i);
    free(k);
    free(l);
    free(y);
}