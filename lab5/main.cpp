#include <iostream>
#include <cmath>
#include <ctime>
#include <pthread.h>
#include <exception>
#include <array>
#include <vector>

#ifdef __unix__
#include <mpi.h>
#elif defined(_WIN32) || defined(WIN32)
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#endif

#define ITERATIONS_NUM (3)

#define TASK_LIST_SIZE (200)
#define TASKS_TO_SEND (TASK_LIST_SIZE / (10 * (rank + 1))

#define TASK_END (11)
#define TASK_NEED (22)
#define TASK_SEND (33)
#define TASK_RECEIVED (44)
#define WORK_DONE (55)

using namespace std;

class Task{
private:
    int repeatNum = 0;
public:
    explicit Task(const int& repeatNum){
        this->repeatNum = repeatNum;
    }
    Task(){ repeatNum = 0; }
    int getRepeatNum() const{ return repeatNum; }
};

class TaskList{
private:
    array<Task, TASK_LIST_SIZE> tasks = {Task(0)};

public:
    TaskList() = default;
    Task& at(const size_t& idx){ return this->tasks.at(idx); }
    void addTask(const Task& task, const int& idx){ tasks.at(idx) = task; }
    size_t getTasksNum(){ return tasks.size(); }
};

void generateTasks(TaskList& taskList, size_t& iterationCnt, pthread_mutex_t& taskListMutex, const int& rank, const int& size){
    pthread_mutex_lock(&taskListMutex);
    for (int i = rank * TASK_LIST_SIZE; i < (rank + 1) * TASK_LIST_SIZE; i++){
        Task newTask(static_cast<int>(fabs(50 - i % TASK_LIST_SIZE) * fabs(rank - (iterationCnt % size)) * 1200));
        taskList.addTask(newTask, i);
    }
    pthread_mutex_unlock(&taskListMutex);
}

void calculate(long double& calculations, TaskList& taskList, pthread_mutex_t& taskListMutex, const size_t& idx){
    pthread_mutex_lock(&taskListMutex);
    for(size_t i = 0; i < taskList.at(idx).getRepeatNum(); ++i){
        calculations += sin(sin(static_cast<double>(i)));
    }
    pthread_mutex_unlock(&taskListMutex);
}

void initMPI(int* rank, int* size){
    MPI_Comm_size(MPI_COMM_WORLD, size);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
}


int main(int argc, char* argv[]) {
    int provided;
    pthread_attr_t attrs;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE){
        std::cout << "Can't init thread" << std::endl;
        MPI_Finalize();
        return 0;
    }
    int size = 0, rank = 0;
    initMPI(&rank, &size);

    double timeStart = 0, timeFinish = 0;
    if(rank == 0) timeStart = MPI_Wtime();

    TaskList taskList;
    pthread_mutex_t taskListMutex;
    long double calculations = 0;
    size_t iterationCnt = 0;




    if(rank == 0) {
        timeFinish = MPI_Wtime();
        cout << "Time: " << timeFinish - timeStart << endl;
    }
    MPI_Finalize();
    return 0;
}
