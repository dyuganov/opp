#include<mpi.h> // Подключение библиотеки MPI
#include<stdio.h>

int main(int argc, char* argv[])
{
	int size, rank;
	MPI_Init(&argc, &argv); // Инициализация MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size); // Получение числа процессов
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Получение номера процесса
	printf("Hello from process # %d of %d\n", rank, size);
	MPI_Finalize(); // Завершение работы MPI
	return 0;
}