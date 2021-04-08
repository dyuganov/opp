#include <iostream>
#include "Matrix.h"

int main(int argc, char *argv[]) {
    int size, rank, part_size, extended; //количество процессов, номер текущего процесса, размер части матрицы, максимальная длина для векторов b и x
    int* lengths; //массив размеров частей матрицы в строках
    int* displs; // массив отступов от 0-й строки матрицы до каждой из частей
    double* partRes; // частичный вектор-результат вычислений в цикле
    double* matrix; // исходная матрица
    double* b;//частичный вектор b
    double* x;//частичный вектор x
    double* part_mat; //отрезанная часть матрицы
    double collectedResult[N], vectorBNorm, sumResult; //вектор для сбора векторов частичных векторов b, длина целого вектора b, число для сбора всех Ax-b с процессов
    double* full_x; //полный вектор x
    double* full_b; //полный вектор b

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //блок инициализации процессов в количестве size в области связи (коммуникаторе) MPI_COMM_WORLD
    //ниже блок выделения памяти и вычисления вышеописанных переменных
    lengths = (int*)malloc(size* sizeof(int));
    getAllLengths(size, lengths);
    displs = (int*)malloc(size*sizeof(int));
    prepareDispls(size, lengths, displs);
    part_size = lengths[rank];
    extended = (int)ceil((double)N / size);
    partRes = (double*)calloc(extended, sizeof(double));
    matrix = (double*)malloc(N * N * sizeof(double));
    b = (double*)malloc(extended * sizeof(double));
    x = (double*)malloc(extended * sizeof(double));
    full_x = (double*)malloc(N * sizeof(double));
    full_b = (double*)malloc(N * sizeof(double));
    part_mat = (double*)calloc(part_size * N , sizeof(double));
    //для MPI_Scatterv 2-м и 3-м аргументами нужны массивы с КОЛИЧЕСТВОМ ЭЛЕМЕНТОВ, поэтому
    //здесь создаются 2 массива под Кол-во эл-ов в строках отсеченной части и кол-вом эл-ов на которые надо сдвинуться до след. части
    int* sendLengths = (int*)malloc(size * sizeof(int));
    int* sendDispls = (int*)malloc(size * sizeof(int));

    if (rank == 0) {
        initBX(full_x, full_b);
        matrixInit(matrix);
        vectorBNorm = getVecLen(full_b, N);
    }
    MPI_Bcast(&vectorBNorm,1, MPI_DOUBLE,0, MPI_COMM_WORLD);

    MPI_Scatterv(full_x, lengths, displs , MPI_DOUBLE, x,
                 extended, MPI_DOUBLE, 0 , MPI_COMM_WORLD); //рассылка частей вектора full_x в x

    MPI_Scatterv(full_b, lengths, displs , MPI_DOUBLE, b,
                 extended, MPI_DOUBLE, 0 , MPI_COMM_WORLD); //рассылка частей вектора b

    double time = MPI_Wtime(); //замер времени для последующего вывода на экран

    for (int i = 0; i < size; i++) { //инициализация массивов длин и сдвигов в элементах
        sendLengths[i] = lengths[i] * N;
        sendDispls[i] = displs[i] * N;
    }

    MPI_Scatterv(matrix, sendLengths, sendDispls , MPI_DOUBLE, part_mat,
                 N * part_size, MPI_DOUBLE, 0 , MPI_COMM_WORLD); //рассылка матрицы matrix на части part_mat

    while (true){
        mulMatrixPart(part_mat, x, partRes, length s, displs, size, rank); // Ax
        vecSub(partRes, b, partRes, extended); //Ax - b
        double sum = 0;
        for (size_t i = 0; i < part_size; i++) { //находим (Ax-b)^2 для текущего процесса
            sum += partRes[i] * partRes[i];
        }
        MPI_Allreduce(&sum, &sumResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // суммируем все (Ax-b)^2 для проверки
        if (isEnd(sumResult, vectorBNorm)) // проверка на ( ||Ax-b|| / ||b|| )< epsilon
            break;
        vecMulScalar(partRes, t, partRes, extended); // t(Ax -b)
        vecSub(x, partRes, x, extended); //x^2 = x - t(Ax - b)
}

    if (rank == 0)
        std::cout << "final time is: " << MPI_Wtime() - time << std::endl;
    printf("pr(%d): ", rank);
    printVector(x,  part_size);
    freeAll(matrix, b, x, part_mat, lengths, sendLengths,
            sendDispls, displs, full_x, full_b, partRes);
    MPI_Finalize();
    return 0;
}
