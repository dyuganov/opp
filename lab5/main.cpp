#include <iostream>
#include <cmath>
#include <pthread.h>

#ifdef __unix__
#include <mpi.h>
#elif defined(_WIN32) || defined(WIN32)
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#endif

#define ITERATIONS_NUM (3)

#define TASK_LIST_SIZE (250)
#define TASKS_TO_SEND (TASK_LIST_SIZE/(10*(procRank+1)))//

#define NO_TASK (11)
#define TASK_NEED (22)
#define TASK_SEND (33)
#define TASK_RECEIVED (44)
#define WORK_DONE (55)

#define REQUEST_TAG (101)
#define ANSWER_TAG (102)

using namespace std;

/*class Task{
private:
    int repeatNum = 0;
public:
    explicit Task(const int& repeatNum){
        this->repeatNum = repeatNum;
    }
    Task(){ repeatNum = 0; }
    int getRepeatNum() const{ return repeatNum; }
};*/

pthread_mutex_t taskListMutex;
long double calculations = 0;
size_t iterationCnt = 0;
int size = 0, procRank = 0;
int ownCurrentTaskIdx, givenTask;
int givenTasksCnt, receivedTasksCnt;
int recvMsg_1, sendMsg_1, recvMsg_2, sendMsg_2;

class TaskList{
private:
    int* tasks = nullptr;
public:
    TaskList() = default;
    ~TaskList() {
        delete[] tasks;
        tasks = nullptr;
    }
    void initList() { tasks = new int[TASK_LIST_SIZE * size]; }
    int& at(const size_t& idx){ return this->tasks[idx]; }
    void addTask(const int& task, const size_t& idx){ if(idx < TASK_LIST_SIZE * size) tasks[idx] = task; }
    int* getArray(){ return tasks; }
};
TaskList taskList;

bool isJoinError(const int& joinResult){
    if(joinResult != 0){
        fprintf(stderr, "pthread_join error\n");
        return true;
    }
    return false;
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

void generateTasks() {
    pthread_mutex_lock(&taskListMutex);
    for (int i = procRank * TASK_LIST_SIZE; i < (procRank + 1) * TASK_LIST_SIZE; i++) {
        taskList.at(i) = fabs(50 - i % TASK_LIST_SIZE) * fabs(procRank - (iterationCnt % size)) * 1200;
    }
    pthread_mutex_unlock(&taskListMutex);
}

void calculate(int idx) {
    for (size_t i = 0; i < taskList.at(idx); i++){
        calculations += sin(sin(i));
    }
}

int getNewTask(const int& provider) {
    sendMsg_1 = TASK_NEED;

    MPI_Send(&sendMsg_1, 1, MPI_INT, provider, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&recvMsg_1, 1, MPI_INT, provider, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (recvMsg_1 == NO_TASK) return false;
    if (recvMsg_1 == TASK_SEND) {
        pthread_mutex_lock(&taskListMutex);
        MPI_Recv(&receivedTasksCnt, 1, MPI_INT, provider, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(taskList.getArray()[provider * TASK_LIST_SIZE]), receivedTasksCnt, MPI_INT, provider, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //MPI_Recv(&(taskList.at(provider * TASK_LIST_SIZE)), receivedTasksCnt, MPI_INT, provider, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        givenTask = provider * TASK_LIST_SIZE;
        pthread_mutex_unlock(&taskListMutex);
        return true;
    }
    sendMsg_1 = TASK_RECEIVED;
    MPI_Send(&sendMsg_1, 1, MPI_INT, provider, REQUEST_TAG, MPI_COMM_WORLD);
}

void* resolveTask(void* args) {
    generateTasks();
    int tasksCounter = 0;
    int totalCounter = 0;

    double avgImbalance = 0;

    while (iterationCnt < ITERATIONS_NUM) {
        if (procRank == 0) std::cout << "__________ ITERATION # " << iterationCnt << " __________" << std::endl;
        double timeStart = MPI_Wtime();
        pthread_mutex_lock(&taskListMutex);
        ownCurrentTaskIdx = TASK_LIST_SIZE * procRank;
        givenTasksCnt = 0;
        pthread_mutex_unlock(&taskListMutex);

        while (ownCurrentTaskIdx < TASK_LIST_SIZE * (procRank + 1) - givenTasksCnt){
            calculate(ownCurrentTaskIdx);
            pthread_mutex_lock(&taskListMutex);
            ownCurrentTaskIdx++;
            pthread_mutex_unlock(&taskListMutex);
            tasksCounter++;
        }

        while (true) {
            for (int pr = 0; pr < size; ++pr){
                if (pr != procRank && getNewTask(pr)){
                    for (int i = 0; i < receivedTasksCnt; i++){
                        calculate(givenTask);
                        tasksCounter++;
                        givenTask++;
                    }
                    break;
                }
            }
            break;
        }

        double timeFinish = MPI_Wtime();
        double maxTime, minTime;
        std::cout << "Time: " << (timeFinish - timeStart) << " sec. Rank: " << procRank << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        double time = (timeFinish - timeStart);
        MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&time, &minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&tasksCounter, &totalCounter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (procRank == 0) {
            std::cout << "Total tasks: " << totalCounter << std::endl;
            std::cout << "Imbalance: " << maxTime - minTime << std::endl;
            std::cout << "Imbalance / MAX: " << (maxTime - minTime) / maxTime * 100 << " %" << std::endl;
            avgImbalance += (maxTime - minTime) / maxTime * 100;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << tasksCounter << " tasks - rank: " << procRank << std::endl;
        std::cout << calculations << " - current result - rank: " << procRank << std::endl;
        generateTasks();
        tasksCounter = 0;
        iterationCnt++;
    }

    pthread_mutex_lock(&taskListMutex);
    sendMsg_1 = WORK_DONE;
    MPI_Send(&sendMsg_1, 1, MPI_INT, procRank, REQUEST_TAG, MPI_COMM_WORLD);
    pthread_mutex_unlock(&taskListMutex);
    MPI_Barrier(MPI_COMM_WORLD);
    if (procRank == 0) std::cout << "Avg imbalance: " << avgImbalance / ITERATIONS_NUM << " %" << std::endl;
    return nullptr;
}

void* sendTask(void* args){
    MPI_Status receiver;
    while (iterationCnt < ITERATIONS_NUM){
        MPI_Recv(&recvMsg_2, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &receiver);
        if (receiver.MPI_SOURCE != procRank && recvMsg_2 == TASK_NEED) {
            if (ownCurrentTaskIdx >= (procRank + 1) * TASK_LIST_SIZE - givenTasksCnt - 1){
                sendMsg_2 = NO_TASK;
                MPI_Send(&sendMsg_2, 1, MPI_INT, receiver.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
                continue;
            }
            int* taskForReceiver;
            int extraTasks = 0;
            pthread_mutex_lock(&taskListMutex);
            givenTasksCnt += TASKS_TO_SEND;
            if (givenTasksCnt + ownCurrentTaskIdx >= (procRank + 1) * TASK_LIST_SIZE - 1){
                extraTasks = givenTasksCnt + ownCurrentTaskIdx - (procRank + 1) * TASK_LIST_SIZE + 1;
                givenTasksCnt = (procRank + 1) * TASK_LIST_SIZE - 1 - ownCurrentTaskIdx;
            }
            taskForReceiver = &taskList.getArray()[(procRank + 1) * TASK_LIST_SIZE - givenTasksCnt - 1];
            pthread_mutex_unlock(&taskListMutex);
            sendMsg_2 = TASK_SEND;
            int tasksToSend = TASKS_TO_SEND - extraTasks;
            MPI_Send(&sendMsg_2, 1, MPI_INT, receiver.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            MPI_Send(&tasksToSend, 1, MPI_INT, receiver.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            MPI_Send(taskForReceiver, tasksToSend, MPI_INT, receiver.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            continue;
        }
        else if (recvMsg_2 == TASK_RECEIVED) continue;
        else if (recvMsg_2 == WORK_DONE) break;
    }
    return nullptr;
}

int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE){
        std::cout << "Can't init thread" << std::endl;
        MPI_Finalize();
        return 0;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    taskList.initList();

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

    pthread_t sendingThread;
    pthread_t resolvingTasksThread;
    pthread_create(&sendingThread, &attrs, sendTask, nullptr);
    pthread_create(&resolvingTasksThread, &attrs, resolveTask, nullptr);
    int joinResult1 = pthread_join(sendingThread, nullptr);
    int joinResult2 = pthread_join(resolvingTasksThread, nullptr);
    if(isJoinError(joinResult1) || isJoinError(joinResult2)) {
        pthread_attr_destroy(&attrs);
        pthread_mutex_destroy(&taskListMutex);
        MPI_Finalize();
    }

    pthread_attr_destroy(&attrs);
    pthread_mutex_destroy(&taskListMutex);
    MPI_Finalize();
    return 0;
}
