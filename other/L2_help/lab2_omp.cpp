#include <cmath>
#include <ctime>
#include <iostream>
#include <memory.h>
#include <omp.h>
#include <sys/time.h>
#include <xmmintrin.h>

#define N 2048
#define iterations 10
#define NUM_THREADS 16
using namespace std;

void printVector(float* A)
{
    for (size_t i = 0; i < N * N; i++) {
        if (!(i % N))
            cout << endl;
        cout << A[i] << " ";
    }
}

void mulDef(const float* A, const float* B, float* BA)
{
    for (int i = 0; i < N; i++) {
        float* ba = BA + i * N;
        for (int j = 0; j < N; j++) {
            ba[j] = 0;
        }
        for (int k = 0; k < N; k++) {
            const float* b = B + k * N;
            float a = A[i * N + k];

            for (int j = 0; j < N; k++) {
                ba[j] += a * b[j];
            }
        }
    }
}

void mulSSE(const float* B, const float* A, float* C)
{
#pragma omp parallel
    {
        float* sum_of_lines = new float[N];
        float* tmp = new float[N];

#pragma omp for private(i, j, sum_of_lines, tmp) shared(A, B, C)
        for (int i = 0; i < N; i++)
        {
            memset(sum_of_lines, 0, N * sizeof(float));
            for (int j = 0; j < N; j++)
            {
                __m128 element_of_B = _mm_set1_ps(B[i * N + j]);

                for (int k = 0; k < N; k += 4)
                {
                    __m128 part_of_A = _mm_load_ps(A + (j * N + k));
                    __m128 mul_of_lines = _mm_mul_ps(element_of_B, part_of_A);
                    _mm_store_ps(tmp + k, mul_of_lines);
                }

                for (int p = 0; p < N; p += 4)
                {
                    __m128 sum_of_lines_part = _mm_load_ps(sum_of_lines + p);
                    __m128 tmp_2 = _mm_load_ps(tmp + p);
                    __m128 tmp_lines_sum = _mm_add_ps(sum_of_lines_part, tmp_2);
                    _mm_store_ps(sum_of_lines + p, tmp_lines_sum);
                }
            }
            memcpy(C + i * N, sum_of_lines, sizeof(float) * N);
        }
        delete[] sum_of_lines;
        delete[] tmp;
    }
}

void sumSSE(float* A, float* B, float* C, bool sign)
{
    __m128 p;

#pragma omp parallel for private(p)
    for (int i = 0; i < N; i++)
    {
        __m128* xx;
        __m128* yy;
        xx = (__m128*)(A + i * N);
        yy = (__m128*)(B + i * N);
        for (int k = 0; k < N / 4; k++)
        {
            if (sign)
                p = _mm_add_ps(xx[k], yy[k]);
            else
                p = _mm_sub_ps(xx[k], yy[k]);
            _mm_store_ps(C + i * N + k * 4, p);
        }
    }
}

void cpy(float* A, float* B)
{
    int i, j;
#pragma omp parallel for private(i, j) shared(A, B)
    for (i = 0; i < N; i++) {
        float* a = A + i * N;
        float* b = B + i * N;
        for (j = 0; j < N; j++)
            a[j] = b[j];
    }
}

double SSEMatrixInverse(float* matrix, float* result)
{

    float* I = new float[N * N];
    float* buff = new float[N * N];
    float* B = new float[N * N];
    float* BA = new float[N * N];
    float* tmp = new float[N * N];
    float* R = new float[N * N];

    int i, j;
    {
        float maxRowSum = 0;
        float maxColumnSum = 0;
        float rowSum = 0;
        float columnSum = 0;
#pragma omp parallel for private(i, j) shared(matrix, rowSum, columnSum, I, buff, maxColumnSum, maxRowSum)
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++)

            {
                rowSum += fabs(matrix[N * i + j]);
                columnSum += fabs(matrix[j * N + i]);
                I[N * i + j] = (i == j);
                buff[N * i + j] = (i == j);
            }
            if (columnSum > maxColumnSum) {
                maxColumnSum = columnSum;
            }
            if (rowSum > maxRowSum) {
                maxRowSum = rowSum;
            }
            rowSum = 0;
            columnSum = 0;
        }
        float max = maxColumnSum * maxRowSum;
#pragma omp parallel for private(i, j) shared(B, matrix, max)
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                B[N * i + j] = matrix[j * N + i] / max;
            }
        }
    }
    mulSSE(B, matrix, BA);
    sumSSE(I, BA, R, false);

    for (i = 0; i < iterations - 1; i++)
{
        mulSSE(I, R, tmp);
        cpy(I, tmp);
        sumSSE(buff, I, buff, true);
    }
    mulSSE(buff, B, result);
    clock_t c_end = clock();

    delete[] I;
    delete[] B;
    delete[] A;
    delete[] tmp;

    delete[] R;
    delete[] buff;
    return 1;
}

int main()
{
    omp_set_num_threads(NUM_THREADS);

    struct timeval start, end;

    float* matrix = new float[N * N];
    float* inversedMatrix = new float[N * N];
    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; ++i) {
            matrix[k * N + i] = rand() % 10;
            inversedMatrix[k * N + i] = 0;
        }
    }

    gettimeofday(&start, NULL);
    double res = SSEMatrixInverse(matrix, inversedMatrix);
    gettimeofday(&end, NULL);

    cout << "SSE Inversion : " << end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec) << "s" << endl;

    float* result = new float[N * N];

    //mulDef(matrix, inversedMatrix, result);

     //printVector(result);

    delete[] matrix;
    delete[] inversedMatrix;
    return 0;
}
