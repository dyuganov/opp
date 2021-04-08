#include "stdio.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstring>
#include "mpi.h"

#define N 3200

#define EPS 0.00001

void printVector(double* a)
{
	std::cout << "[ ";

	for (size_t i = 0; i < N; i++)
	{
		std::cout << a[i] << " ";
	}

	std::cout << "]" << std::endl;
}

void printVector(int* a)
{
	std::cout << "[ ";

	for (size_t i = 0; i < N; i++)
	{
		std::cout << a[i] << " ";
	}

	std::cout << "]" << std::endl;
}

void printMatrix(double* a)
{
	std::cout << std::endl;

	std::cout << "[ ";

	for (size_t i = 0; i < N * N; i++)
	{
		std::cout << a[i] << " ";

		if (i % N == 0)
		{
			std::cout << std::endl;
		}
	}

	std::cout << "]" << std::endl;
	std::cout << std::endl;
}



void sumVector(double* a, double* b, double* c, bool sign, int vector_size)
{
	if (sign)
	{
		for (size_t i = 0; i < vector_size; i++)
		{
			c[i] = a[i] + b[i];
		}
	}
	else
	{
		for (size_t i = 0; i < vector_size; i++)
		{
			c[i] = a[i] - b[i];
		}
	}
}

double dotProduct(double* a, double* b, int vector_size)
{
	double dp = 0;

	for (size_t i = 0; i < vector_size; i++)
	{
		dp += a[i] * b[i];
	}

	return dp;
}

void mulScalar(double scalar, double* a, double* b, int vector_size)
{
	for (size_t i = 0; i < vector_size; i++)
	{
		b[i] = scalar * a[i];
	}
}

void mulMatrixAndVector(double* matrix, double* a, double* b)
{
	for (size_t i = 0; i < N; i++)
	{
		b[i] = 0;

		for (size_t j = 0; j < N; j++)
		{
			b[i] += a[j] * matrix[i * N + j];
		}
	}
}

int mulRowAndVector(double* row, double* vector, double* resPart, int rank, int* block)
{
	int d = 0, q = 0;
	for (int i = 0; i < block[rank]; i++) {
		resPart[q] += row[i] * vector[d];
		d++;
		if (d == N)
		{
			q++;
			d = 0;
		}
	}
	return q;

}

void mulColumnAndVector(double* column, double* vector, double* resPart, int rank, int* block)
{
	int d = -1;
	for (int i = 0; i < block[rank]; i++) 
	{
		if (i % N == 0)
		{
			d++;
		}

		resPart[i % N] += column[i] * vector[d];

	}
}

double getVectorLength(double* a)
{
	double length = 0;

	for (size_t i = 0; i < N; i++)
	{
		length += a[i] * a[i];
	}

	length = sqrt(length);

	return length;
}

double getVectorPartLength(double* a, int* vector_parts_capacity, int rank)
{
	double length_sqr = 0;

	for (size_t i = 0; i < vector_parts_capacity[rank]; i++)
	{
		length_sqr += a[i] * a[i];
	}

	return length_sqr;
}

void generateMatrix(double* matrix)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = i; j < N; j++)
		{
			matrix[i * N + j] = rand() % 10;
			matrix[j * N + i] = matrix[i * N + j];
		}
	}

	for (size_t i = 0; i < N; i++)
	{
		matrix[i * N + i] += 300;
	}
}

void generateVector(double* vector)
{
	for (size_t i = 0; i < N; i++)
	{
		vector[i] = rand() % 10;
	}
}

void initArrays(int procNum, int second_N, int* matrix_parts_capacity, int* matrix_index_offset, int* vector_parts_capacity, int* vector_index_offset)
{
	int rem = second_N % procNum;
	second_N -= rem;
	int count = second_N / procNum;

	for (int i = 0; i < procNum; i++)
	{
		matrix_parts_capacity[i] = count * N;
	}
	for (int i = 0; i < rem; i++)
	{
		matrix_parts_capacity[i] += N;
	}

	int sm = 0;
	for (int i = 0; i < procNum; i++)
	{
		matrix_index_offset[i] = sm;
		sm += matrix_parts_capacity[i];
	}

	for (int i = 0; i < procNum; i++)
	{
		vector_parts_capacity[i] = matrix_parts_capacity[i] / N;
	}

	count = 0;
	for (int i = 0; i < procNum; i++)
	{
		vector_index_offset[i] = count;

		count += vector_parts_capacity[i];
	}
}



double* runAlgorythmPARALLEL_2(double* A, double* x, double* b, int rank, int procNum)
{
	// VALUES FOR ALGORYTHM

	b = new double[N]();
	x = new double[N]();
	A = new double[N * N]();
	double* r = new double[N]();
	double* z = new double[N]();
	double b_length = 0;
	double r_length = 0;

	double dp_1 = -1, dp_2;
	double number_of_iterations = 0;

	double alpha = 0;
	double beta = 0;

	double* Ax = new double[N]();

	// VALUES FOR PARALLEL MULTIPLICATION

	int* matrix_parts_capacity = new int[procNum]();
	int* matrix_index_offset = new int[procNum]();
	int* vector_parts_capacity = new int[procNum]();
	int* vector_index_offset = new int[procNum]();
	int second_N = N;

	if (rank == 0)
	{
		initArrays(procNum, second_N, matrix_parts_capacity, matrix_index_offset, vector_parts_capacity, vector_index_offset);

		generateVector(b);
		generateVector(x);
		generateMatrix(A);
	}

	int part_size = (second_N / procNum + 1) * N;
	double* part_of_matrix = new double[part_size]();

	MPI_Bcast(matrix_parts_capacity, procNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(matrix_index_offset, procNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vector_parts_capacity, procNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vector_index_offset, procNum, MPI_INT, 0, MPI_COMM_WORLD);

	int part_vector_size = matrix_parts_capacity[rank] / N;

	double* part_of_vector_1 = new double[N]();
	double* part_of_vector_2 = new double[part_vector_size]();
	double* part_of_vector_3 = new double[part_vector_size]();
	double* part_of_vector_4 = new double[part_vector_size]();
	double* part_of_Az = new double[part_vector_size]();
	double* part_of_alphaz = new double[part_vector_size]();
	double* part_of_betaz = new double[part_vector_size]();

	MPI_Scatterv(x, vector_parts_capacity, vector_index_offset, MPI_DOUBLE, part_of_vector_2, part_vector_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(b, vector_parts_capacity, vector_index_offset, MPI_DOUBLE, part_of_vector_4, part_vector_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(A, matrix_parts_capacity, matrix_index_offset, MPI_DOUBLE, part_of_matrix, part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	mulColumnAndVector(part_of_matrix, part_of_vector_2, part_of_vector_1, rank, matrix_parts_capacity);
	MPI_Reduce_scatter(part_of_vector_1, part_of_vector_3, vector_parts_capacity, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	double sqr_part_length = getVectorPartLength(part_of_vector_4, vector_parts_capacity, rank);
	MPI_Allreduce(&sqr_part_length, &b_length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	b_length = sqrt(b_length);

	sumVector(part_of_vector_4, part_of_vector_3, part_of_vector_4, 0, vector_parts_capacity[rank]);

	MPI_Scatterv(z, vector_parts_capacity, vector_index_offset, MPI_DOUBLE, part_of_vector_3, part_vector_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	memcpy(part_of_vector_3, part_of_vector_4, sizeof(double) * part_vector_size);

	sqr_part_length = getVectorPartLength(part_of_vector_4, vector_parts_capacity, rank);
	MPI_Allreduce(&sqr_part_length, &r_length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	r_length = sqrt(r_length);

	double tmp_eps = 0;
	int counter = 0;

	// z in part_of_vector_3
	// x in part_of_vector_2
	// r in part_of_vector_4

	while (r_length / b_length >= EPS)
	{
		if (r_length / b_length > tmp_eps)
		{
			tmp_eps = r_length / b_length;
			counter++;
		}
		if (counter == 5)
		{
			std::cout << "Can't resolve matrix." << std::endl;
			break;
		}

		number_of_iterations++;

		for (size_t i = 0; i < N; i++)
		{
			part_of_vector_1[i] = 0;
		}

		mulColumnAndVector(part_of_matrix, part_of_vector_3, part_of_vector_1, rank, matrix_parts_capacity);
		MPI_Reduce_scatter(part_of_vector_1, part_of_Az, vector_parts_capacity, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		if (dp_1 == -1)
		{
			double tmp = dotProduct(part_of_vector_4, part_of_vector_4, vector_parts_capacity[rank]);
			MPI_Allreduce(&tmp, &dp_1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		}
		else
		{
			dp_1 = dp_2;
		}

		double tmp_1 = dotProduct(part_of_Az, part_of_vector_3, vector_parts_capacity[rank]);
		MPI_Allreduce(&tmp_1, &dp_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		alpha = dp_1 / dp_2;

		mulScalar(alpha, part_of_vector_3, part_of_alphaz, vector_parts_capacity[rank]);

		sumVector(part_of_vector_2, part_of_alphaz, part_of_vector_2, 1, vector_parts_capacity[rank]);

		mulScalar(alpha, part_of_Az, part_of_Az, vector_parts_capacity[rank]);

		sumVector(part_of_vector_4, part_of_Az, part_of_vector_4, 0, vector_parts_capacity[rank]);

		double tmp_2 = dotProduct(part_of_vector_4, part_of_vector_4, vector_parts_capacity[rank]);
		MPI_Allreduce(&tmp_2, &dp_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		beta = dp_2 / dp_1;

		mulScalar(beta, part_of_vector_3, part_of_betaz, vector_parts_capacity[rank]);
		sumVector(part_of_vector_4, part_of_betaz, part_of_vector_3, 1, vector_parts_capacity[rank]);

		sqr_part_length = getVectorPartLength(part_of_vector_4, vector_parts_capacity, rank);
		MPI_Allreduce(&sqr_part_length, &r_length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		r_length = sqrt(r_length);
	}

	MPI_Gatherv(part_of_vector_2, vector_parts_capacity[rank], MPI_DOUBLE, x, vector_parts_capacity, vector_index_offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	delete[] r;
	delete[] z;
	delete[] part_of_Az;
	delete[] Ax;
	delete[] part_of_alphaz;
	delete[] part_of_betaz;

	delete[] part_of_matrix;
	delete[] matrix_parts_capacity;
	delete[] vector_index_offset;
	delete[] matrix_index_offset;
	delete[] vector_parts_capacity;
	delete[] part_of_vector_4;
	delete[] part_of_vector_1;
	delete[] part_of_vector_2;
	delete[] part_of_vector_3;

	//if (rank == 0)
	//{
	//	double* tmp = new double[N]();
	//	mulMatrixAndVector(A, x, tmp);
	//	printVector(tmp);
	//	printVector(b);
	//	delete[] tmp;
	//}

	if (rank == 0) std::cout << "Number of iterations: " << number_of_iterations << std::endl;
	return x;
}

int main(int argc, char** argv)
{
	int procNum, rank;
	clock_t start, finish;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::cout << "My rank is " << rank << std::endl;

	double* A = nullptr;
	double* b = nullptr;
	double* x = nullptr;

	if (rank == 0) start = clock();

	x = runAlgorythmPARALLEL_2(A, x, b, rank, procNum);

	if (rank == 0) finish = clock();

	if (rank == 0)
	{
		std::cout << "Time: " << (double)((finish - start) / CLOCKS_PER_SEC) << std::endl;
	}

	/*if (rank == 0)
	{
		double* tmp = new double[N]();
		mulMatrixAndVector(A, x, tmp);
		printVector(tmp);
		printVector(b);
		delete[] tmp;
	}*/

	delete[] x;
	delete[] b;
	delete[] A;

	MPI_Finalize();
	return 0;
}
