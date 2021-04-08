#include <cmath>
#include <ctime>
#include <iostream>
#include <xmmintrin.h>
#define N 2048
#define iterations 10
using namespace std;

void mulSSE(const float* B, const float* A, float* C)
{
        float* sum_of_lines = new float[N];
        float* tmp = new float[N];


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
void mulDef(const float* A, const float* B, float* BA)
{
	for (int i = 0; i < N; i++)
	{
		float* ba = BA + i * N;
		for (int j = 0; j < N; j++)
		{
			ba[j] = 0;
		}
		for (int k = 0; k < N; k++)
		{
			const float* b = B + k * N;
			float a = A[i * N + k];

			for (int j = 0; j < N; j++)
			{
				ba[j] += a * b[j];
			}
		}
	}
}
void sumDef(float* A, float* B, float* C, bool sign)
{
	if (sign)
		for (int i = 0; i < N; i++)
		{
			float* a = A + i * N;
			float* b = B + i * N;
			float* c = C + i * N;
			for (int j = 0; j < N; j++)
				c[j] = a[j] + b[j];

		}
	else
		for (int i = 0; i < N; i++)
		{
			float* a = A + i * N;
			float* b = B + i * N;
			float* c = C + i * N;
			for (int j = 0; j < N; j++)
				c[j] = a[j] - b[j];

		}
}
void sumSSE(float* A, float* B, float* C, bool sign)
{
	__m128 p;
	__m128* xx;
	__m128* yy;
	for (int i = 0; i < N; i++)
	{
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

	for (int i = 0; i < N; i++)
	{
		float* a = A + i * N;
		float* b = B + i * N;
		for (int j = 0; j < N; j++)
			a[j] = b[j];

	}
}
double defaultMatrixInverse(float* matrix, float* result)
{
	clock_t c_start = clock();
	float* I = new float[N * N];
	float* buff = new float[N * N];
	float* B = new float[N * N];
	float* BA = new float[N * N];
	float* tmp = new float[N * N];
	float* R = new float[N * N];
	{
		float maxRowSum = 0;
		float maxColumnSum = 0;
		float rowSum = 0;
		float columnSum = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				rowSum += fabs(matrix[N * i + j]);
				columnSum += fabs(matrix[j * N + i]);
				I[N * i + j] = (i == j);
				buff[N * i + j] = (i == j);
			}
			if (columnSum > maxColumnSum)
			{
				maxColumnSum = columnSum;
			}
			if (rowSum > maxRowSum)
			{
				maxRowSum = rowSum;
			}
			rowSum = 0;
			columnSum = 0;
		}
		float max = maxColumnSum * maxRowSum;
		for (int i = 0; i < N; i++)

		{
			for (int j = 0; j < N; j++)
			{
				B[N * i + j] = matrix[j * N + i] / max;
			}
		}
	}
	mulDef(B, matrix, BA);
	sumDef(I, BA, R, false);
	for (size_t i = 0; i < iterations - 1; i++)
	{
		mulDef(I, R, tmp);
		cpy(I, tmp);
		sumDef(buff, I, buff, true);
	}
	mulDef(buff, B, result);
	clock_t c_end = clock();
	double totalTime = 1000.0 * ((double)(c_end - c_start) / CLOCKS_PER_SEC);
	delete[](I);
	delete[](B);
	delete[](BA);
	delete[](tmp);
	delete[](R);
	delete[](buff);
	return totalTime;
}
double SSEMatrixInverse(float* matrix, float* result)
{
	clock_t c_start = clock();
	float* I = new float[N * N];
	float* buff = new float[N * N];
	float* B = new float[N * N];
	float* BA = new float[N * N];
	float* tmp = new float[N * N];
	float* R = new float[N * N];
	{
		float maxRowSum = 0;
		float maxColumnSum = 0;
		float rowSum = 0;
		float columnSum = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)

			{
				rowSum += fabs(matrix[N * i + j]);
				columnSum += fabs(matrix[j * N + i]);
				I[N * i + j] = (i == j);
				buff[N * i + j] = (i == j);
			}
			if (columnSum> maxColumnSum)
			{
				maxColumnSum = columnSum;
			}
			if (rowSum> maxRowSum)
			{
				maxRowSum = rowSum;
			}
			rowSum = 0;
			columnSum = 0;
		}
		float max = maxColumnSum * maxRowSum;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				B[N * i + j] = matrix[j * N + i] / max;
			}
		}
	}
	mulSSE(B, matrix, BA);
	sumSSE(I, BA, R, false);

	for (size_t i = 0; i < iterations - 1; i++)
	{
		mulSSE(I, R, tmp);
		cpy(I, tmp);
		sumSSE(buff, I, buff, true);
	}
	mulSSE(buff, B, result);
	clock_t c_end = clock();
	double totalTime = 1000.0 * ((double)(c_end - c_start) / CLOCKS_PER_SEC);
	delete[](I);
	delete[](B);
	delete[](BA);
	delete[](tmp);

	delete[](R);
	delete[](buff);
	return totalTime;
}

int main()
{
	float* matrix = new float[N * N];
	float* inversedMatrix = new float[N * N];
	for (int k = 0; k < N; ++k)
	{
		for (int i = 0; i < N; ++i)
		{
			matrix[k * N + i] = rand() % 10;
			inversedMatrix[k * N + i] = 0;
		}
	}

	double res = defaultMatrixInverse(matrix, inversedMatrix);
	cout<<"Default Inversion : " << res<< "ms" << endl;

	res = SSEMatrixInverse(matrix, inversedMatrix);
	cout<< "SSE Inversion : " << res<< "ms" << endl;

	delete[](matrix);
	delete[](inversedMatrix);
	return 0;
}
