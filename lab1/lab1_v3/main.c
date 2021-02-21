#include<mpi.h> // ����������� ���������� MPI
#include<stdio.h>

int main(int argc, char* argv[])
{
	int size, rank;
	MPI_Init(&argc, &argv); // ������������� MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size); // ��������� ����� ���������
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ��������� ������ ��������
	printf("Hello from process # %d of %d\n", rank, size);
	MPI_Finalize(); // ���������� ������ MPI
	return 0;
}