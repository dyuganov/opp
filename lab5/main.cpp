#include <iostream>
#include <vector>
#include "pthread.h"

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
    size_t repeatNum = 0;

public:
    explicit Task(const size_t& repeatNum){
        this->repeatNum = repeatNum;
    }
    Task() = delete;
    size_t getRepeatNum() const{ return repeatNum; }
};

class TaskList{
private:
    vector<Task> tasks;

public:
    TaskList() { tasks.resize(TASK_LIST_SIZE); }
    Task at(const size_t& idx){ return this->tasks.at(idx); }
    void createTask(const size_t& repeatNum){ tasks.emplace_back(repeatNum); }
    void addTask(const Task& task){ tasks.push_back(task); }
    size_t getTasksNum(){ return tasks.size(); }
};

int main(int argc, char* argv[]) {
    int size = 0, rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cout << "Hello! My rank: " << rank << " of " << size << endl;

    MPI_Finalize();
    return 0;
}
