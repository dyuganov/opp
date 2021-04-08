#include "Matrix.h"

int main(int argc, char *argv[]) {
    int size, rank; // количество процессов, номер текущего процесса
    int* lengths;//массив длин (в строках) частей разрезанной матрицы
    int* displs;//массив отступов (в строках) до частей разрезанной матрицы
    double* multResult; // для Ax на каждом процессе
    double* matrix; // матрица A
    double* b; // вектор b
    double* x; // вектор x
    double* collectedResult; // промежуточный результат
    double vectorBNorm; // длина вектора b

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);//блок инициализации

    lengths = getAllLengths(size);
    displs = prepareDispls(size, lengths);
    int part_size = lengths[rank]; // длина (в строках) части матрицы на текущем процессе
    auto matrix_part = (double*)malloc(part_size * N * sizeof(double)); // массив под отрезанную часть матрицы
    //блок выделения памяти уже описанных переменных
    matrix = (double*)malloc(N * N * sizeof(double));
    b = (double*)malloc(N * sizeof(double));
    x = (double*)malloc(N * sizeof(double));
    collectedResult = (double*)malloc(N * sizeof(double));
    multResult = (double*)malloc(part_size * sizeof(double));

    if (rank == 0 )//инициализируем матрицу только на 0 процессе
        matrixInit(matrix);
    otherInit(b, x);//векторы инициализируются на каждом процессе
    //для MPI_Scatterv нужны lengths и displs 2-м и 3-м аргументом в ЭЛЕМЕНТАХ, поэтмоу создаем новые массивы для удобства
    auto send_lengths = (int*)malloc(size * sizeof(double));
    auto send_displs = (int*)malloc(size * sizeof(double));

    for(int i = 0; i < size; i++){ //инициализируем их
        send_displs[i] = N * displs[i];
        send_lengths[i] = N * lengths[i];
    }

    vectorBNorm = getVecLen(b);//находим длину вектора b
    double time = MPI_Wtime();//замер времени для дальнейшего вывода на экран
    MPI_Scatterv(matrix, send_lengths, send_displs , MPI_DOUBLE, matrix_part,
                 N * part_size, MPI_DOUBLE, 0 , MPI_COMM_WORLD);//разрезаем матрицу А на части

    while (1) {
        mulMatrixPart(matrix_part, lengths[rank], x, multResult);//Ax
        MPI_Allgatherv(multResult, lengths[rank], MPI_DOUBLE, collectedResult,
                       lengths, displs, MPI_DOUBLE, MPI_COMM_WORLD);//собираем все результаты со всех процессов (Ax) в один вектор collectedResult
        vecSub(collectedResult, b, collectedResult);//Ax - b
        double tmp = getVecLen(collectedResult); // ||Ax - b||
        if (isEnd(tmp, vectorBNorm))//Проверка на ( ||Ax - b|| / ||b|| ) < epsilon
            break;
        vecMulScalar(t, collectedResult, collectedResult); // t(Ax - b)
        vecSub(x, collectedResult, x);// итерируемся по x дальше
    }

    if (rank == 0) //выводим время только на 0 процессе для чистоты консоли
        std::cout << "time passed: " << MPI_Wtime() - time << std::endl;
    //вывод вектора х на каждом из процессов
    printf("pr(%d): ", rank);
    printVector(x);

    freeAll(lengths, displs, multResult, matrix, x, b,
            collectedResult, matrix_part, send_lengths, send_displs);
    MPI_Finalize();
    return 0;
}