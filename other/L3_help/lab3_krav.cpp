#include "stdio.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstring>
#include "mpi.h"

#define N1 2000
#define N2 2000
#define N3 2000
#define XSIZE 8
#define YSIZE 2
#define X 0
#define Y 1

void mulMatrix(double* A, double* B, double* AB, int n1, int n2, int n3)
{
	for (int i = 0; i < n1; i++)
	{
		double* ba = AB + i * n3;

		for (int j = 0; j < n3; j++)
		{
			ba[j] = 0;
		}
		for (int k = 0; k < n2; k++)
		{
			const double* b = B + k * n3;
			float a = A[i * n2 + k];

			for (int j = 0; j < n3; j++)
			{
				ba[j] += a * b[j];
			}
		}
	}
}


void printMatrix(double* a, int n1, int n2)
{
	std::cout << std::endl;

	for (size_t i = 0; i < n1 * n2; i++)
	{
		if (i % n2 == 0)
		{
			std::cout << std::endl;
		}

		std::cout << a[i] << " ";


	}
	std::cout << std::endl;
}


void generateMatrix(double* matrix, int n1, int n2)
{
	for (size_t i = 0; i < n1; i++)
	{
		for (size_t j = 0; j < n2; j++)
		{
			matrix[i * n2 + j] = rand() % 10;
		}
	}
}

void copyMatrix(double* A, double* B, int n1, int n2, int flag)
{
	if (flag)
	{
		for (int i = 0; i < n1; i++)
		{
			for (int j = 0; j < n2 / XSIZE; j++)
			{
				A[i * n2 / XSIZE + j] = B[i * n2 + j];
			}
		}
	}
	else
	{
		for (int i = 0; i < n1 / YSIZE; i++)
		{
			for (int j = 0; j < n2 / XSIZE; j++)
			{
				A[i * n2 + j] = B[i * n2/XSIZE + j];
			}
		}
	}
}

int main(int argc, char** argv)
{
	int proc_num, rank;
	int dims[2] = {XSIZE, YSIZE};
	int periods[2] = { 0, 0 };
	int cords[2] = { 0 ,0 };
	int remain_dims[2] = { 0, 0 };
	int mapping[XSIZE * YSIZE] = { 0 };
	double start, finish;
	MPI_Comm old_comm, new_comm, row_comm, col_comm;
	old_comm = MPI_COMM_WORLD;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(old_comm, &proc_num);
	MPI_Comm_rank(old_comm, &rank);

	if (proc_num != XSIZE * YSIZE)
	{
		std::cout << "procNum != XSIZE * YSIZE" << std::endl;
		MPI_Finalize();
		return 0;
	}

	double* A = nullptr;
	double* B = nullptr;
	double* AB = nullptr;
	double* part_of_A = nullptr;
	double* part_of_B = nullptr;
	double* part_of_AB = nullptr;

	part_of_A = new double[(N1 * N2) / YSIZE]();
	part_of_B = new double[(N2 * N3) / XSIZE]();
	part_of_AB = new double[(N1 * N3) / XSIZE*YSIZE]();

	if (rank == 0)
	{
		A = new double[N1 * N2]();
		B = new double[N2 * N3]();
		AB = new double[N1 * N3]();

		generateMatrix(A, N1, N2);
		generateMatrix(B, N2, N3);

		//printMatrix(A, N1, N2);
		//printMatrix(B, N2, N3);
	}

	if (rank == 0) start = MPI_Wtime();

	MPI_Dims_create(proc_num, 2, dims);
	MPI_Cart_create(old_comm, 2, dims, periods, 1, &new_comm);
	MPI_Comm_rank(new_comm, &rank);
	MPI_Cart_coords(new_comm, rank, 2, cords);

	if (rank != 0)
	{
		MPI_Send(cords, 2, MPI_INT, 0, 30, new_comm);
	}
	if (rank == 0)
	{
		for (int i = 1; i < proc_num; i++)
		{
			MPI_Recv(cords, 2, MPI_INT, i, 30, new_comm, MPI_STATUS_IGNORE);
			mapping[i] = cords[X] * 10 + cords[Y];
		}

		cords[X] = 0;
		cords[Y] = 0;
	}

	remain_dims[Y] = 0;
	remain_dims[X] = 1; // dims = {1, 0}
	MPI_Cart_sub(new_comm, remain_dims,  &row_comm);

	remain_dims[Y] = 1;
	remain_dims[X] = 0; // dims = {0, 1}
	MPI_Cart_sub(new_comm, remain_dims, &col_comm);

	if (cords[X] == 0)
	{
		MPI_Scatter(A, (N1*N2)/ YSIZE, MPI_DOUBLE, part_of_A, (N1 * N2) / YSIZE, MPI_DOUBLE, 0, col_comm);
	}

	MPI_Bcast(part_of_A, (N1 * N2) / YSIZE, MPI_DOUBLE, 0, row_comm);

	MPI_Datatype COLUMN;
	MPI_Type_vector(N2, N3 / XSIZE, N3, MPI_DOUBLE, &COLUMN);
	MPI_Type_commit(&COLUMN);

	if (rank == 0)
	{
		copyMatrix(part_of_B, B, N2, N3, 1);

		for (int i = 1; i < XSIZE; i++)
		{
			MPI_Send(B + i * (N3 / XSIZE), 1, COLUMN, i, 10, row_comm);
		}
	}

	if (cords[Y] == 0 && rank != 0)
	{
		MPI_Recv(part_of_B, N2 * (N3 / XSIZE), MPI_DOUBLE, 0, 10, row_comm, MPI_STATUS_IGNORE);
	}

	MPI_Bcast(part_of_B, (N2 * N3) / XSIZE, MPI_DOUBLE, 0, col_comm);



	mulMatrix(part_of_A, part_of_B, part_of_AB, N1/YSIZE, N2, N3/XSIZE);

	MPI_Datatype MINOR;
	MPI_Type_vector(N1 / YSIZE, N3 / XSIZE, N3, MPI_DOUBLE, &MINOR);
	MPI_Type_commit(&MINOR);

	if (rank != 0)
	{
		MPI_Send(part_of_AB, (N1 * N3) / (XSIZE * YSIZE), MPI_DOUBLE, 0, 20, new_comm);
	}

	if (rank == 0)
	{
		for (int i = 1; i < proc_num; i++)
		{
			cords[X] = mapping[i] / 10;
			cords[Y] = mapping[i] - cords[X] * 10;
			MPI_Recv(AB + cords[X] * (N3 / XSIZE) + cords[Y] * (N3 * (N1 / YSIZE)), 1, MINOR, i, 20, new_comm, MPI_STATUS_IGNORE);
		}

		copyMatrix(AB, part_of_AB, N1, N3, 0);

		//printMatrix(AB, N1, N3);
	}

	/*MPI_Barrier(old_comm);
	if (rank == 0)
	{

		mulMatrix(A, B, AB, N1, N2, N3);
		printMatrix(AB, N1, N3);
	}*/

	if (rank == 0) finish = MPI_Wtime();

	if (rank == 0)
	{
		std::cout << "Time: " << (finish - start) << std::endl;
	}

	delete[] B;
	delete[] AB;
	delete[] A;
	delete[] part_of_A;
	delete[] part_of_B;
	delete[] part_of_AB;

	MPI_Type_free(&COLUMN);
	MPI_Type_free(&MINOR);
	MPI_Finalize();
	return 0;
}
