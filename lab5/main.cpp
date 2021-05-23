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

#define REQUEST_TAG (101)
#define ANSWER_TAG (102)

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

TaskList taskList;
pthread_mutex_t taskListMutex;
long double calculations = 0;
size_t iterationCnt = 0;

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

bool isAttrError(const int& attrResult){
    if(attrResult != 0){
        fprintf(stderr, "Can't init attrs\n");
        MPI_Finalize();
        return true;
    }
    return false;
}

bool isMutexError(const int& mutexResult){
    if(mutexResult != 0){
        fprintf(stderr, "Can't init mutex\n");
        MPI_Finalize();
        return true;
    }
    return false;
}

int getNewTask(int supplier){
    int msgToSend = TASK_NEED;
    int msgToReceive = 0;

    MPI_Send(&msgToSend, 1, MPI_INT, supplier, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&msgToReceive, 1, MPI_INT, supplier, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if(msgToReceive == TASK_END) return 0;

    if(msgToReceive == TASK_SEND){
        int receiveTasksCnt = 0;
        pthread_mutex_lock(&taskListMutex);
        MPI_Recv(&receiveTasksCnt, 1, MPI_INT, supplier, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&tasklist[sponsor * TASKS_IN_LIST], receiveTasksCnt, MPI_INT, supplier, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        given_task = supplier * TASKS_IN_LIST; // разобраться в логике передачи, исправить
        pthread_mutex_unlock(&taskListMutex);
        return 1;
    }
    msgToSend = TASK_RECEIVED;
    MPI_Send(&msgToSend, 1, MPI_INT, supplier, REQUEST_TAG, MPI_COMM_WORLD);
}

void* resolveTask(void* args){

}

void* sendTask(void* args){

}

int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE){
        std::cout << "Can't init thread" << std::endl;
        MPI_Finalize();
        return 0;
    }
    int size = 0, rank = 0;
    initMPI(&rank, &size);

    pthread_attr_t attrs;
    if(isMutexError(pthread_mutex_init(&taskListMutex, nullptr))) {
        pthread_attr_destroy(&attrs);
        pthread_mutex_destroy(&taskListMutex);
        return 0;
    }
    if(isAttrError(pthread_attr_init(&attrs))) {
        pthread_attr_destroy(&attrs);
        pthread_mutex_destroy(&taskListMutex);
        return 0;
    }
    if(isAttrError(pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE))) {
        pthread_attr_destroy(&attrs);
        pthread_mutex_destroy(&taskListMutex);
        return 0;
    }

    double timeStart = 0, timeFinish = 0;
    if(rank == 0) timeStart = MPI_Wtime();

    pthread_t sendingThread;
    pthread_t resolvingTasksThread;
    pthread_create(&sendingThread, &attrs, sendTask, nullptr);
    pthread_create(&resolvingTasksThread, &attrs, resolveTask, nullptr);
    pthread_join(sendingThread, nullptr);
    pthread_join(resolvingTasksThread, nullptr);

    if(rank == 0) {
        timeFinish = MPI_Wtime();
        cout << "Time: " << timeFinish - timeStart << endl;
    }
    pthread_attr_destroy(&attrs);
    pthread_mutex_destroy(&taskListMutex);
    MPI_Finalize();
    return 0;
}
