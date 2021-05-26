#include "stdio.h"
#include <iostream>
#include <cmath>
#include <ctime>
#define _TIMESPEC_DEFINED
#include <pthread.h>
#include <cstring>
#include "mpi.h"

#define TASKS_IN_LIST 250

#define REQUEST_TAG 100
#define ANSWER_TAG 200

#define NEED_NEW_TASK 111
#define NO_TASK 222
#define SEND_TASK 333
#define TASK_RECEIVED 444
#define WORK_DONE 666

#define WORKING_ITERATIONS 3
#define MUTEX_LOCK pthread_mutex_lock(&mutex)
#define MUTEX_UNLOCK pthread_mutex_unlock(&mutex)
#define TASKS_TO_GIVE TASKS_IN_LIST/(10*(rank+1))

pthread_mutex_t mutex;
int count_of_received_tasks, rank, proc_num, t1_message_to_send, t1_message_to_recv, t2_message_to_send, t2_message_to_recv;
double var_for_calculations;
int iteration_counter = 0;
double start, finish;
int* tasklist;
int current_own_task, given_task;
int count_of_given_tasks;

void calculate(int index)
{
    for (size_t i = 0; i < tasklist[index]; i++)
    {
        var_for_calculations += sin(sin(i));
    }
}

void update_task()
{
    MUTEX_LOCK;// тасклист на всех процесса один большой, но данные хранятся только в нужном куске (плохо!)
    for (int i = rank * TASKS_IN_LIST; i < (rank + 1) * TASKS_IN_LIST; i++)
    {
        tasklist[i] = fabs(50 - i % TASKS_IN_LIST) * fabs(rank - (iteration_counter % proc_num)) * 1200;
    }
    MUTEX_UNLOCK;
}

int get_new_task(int sponsor)
{
    t1_message_to_send = NEED_NEW_TASK;

    MPI_Send(&t1_message_to_send, 1, MPI_INT, sponsor, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&t1_message_to_recv, 1, MPI_INT, sponsor, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (t1_message_to_recv == NO_TASK)
    {
        return 0;
    }
    if (t1_message_to_recv == SEND_TASK)
    {
        MUTEX_LOCK;

        MPI_Recv(&count_of_received_tasks, 1, MPI_INT, sponsor, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&tasklist[sponsor * TASKS_IN_LIST], count_of_received_tasks, MPI_INT, sponsor, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        given_task = sponsor * TASKS_IN_LIST;
        MUTEX_UNLOCK;
        return 1;
    }
    t1_message_to_send = TASK_RECEIVED;
    MPI_Send(&t1_message_to_send, 1, MPI_INT, sponsor, REQUEST_TAG, MPI_COMM_WORLD);
}

void* resolve_task(void* args)
{
    update_task();
    int tasks_counter = 0;
    int total_counter = 0;
    double time = 0;
    double max_time = 0;
    double min_time = 0;
    double avarege_disb = 0;

    while (iteration_counter < WORKING_ITERATIONS)
    {
        if (rank == 0) std::cout << "__________ ITERATION # " << iteration_counter << " __________" << std::endl;
        start = MPI_Wtime();
        MUTEX_LOCK;
        current_own_task = TASKS_IN_LIST * rank;
        count_of_given_tasks = 0;
        MUTEX_UNLOCK;

        while (current_own_task < TASKS_IN_LIST * (rank + 1) - count_of_given_tasks)
        {
            calculate(current_own_task);
            MUTEX_LOCK;
            current_own_task++;
            MUTEX_UNLOCK;
            tasks_counter++;
        }

        bool flag = true;
        while (flag)
        {
            flag = false;
            for (int m = 0; m < proc_num; m++)
            {
                if (m != rank && get_new_task(m))
                {
                    for (int i = 0; i < count_of_received_tasks; i++)
                    {
                        calculate(given_task);
                        tasks_counter++;
                        given_task++;
                    }
                    MUTEX_LOCK;
                    MUTEX_UNLOCK;
                    flag = true;
                    break;
                }
            }
        }

        finish = MPI_Wtime();
        std::cout << "Time: " << (finish - start) << " sec. Rank: " << rank << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        time = (finish - start);
        MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&tasks_counter, &total_counter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            std::cout << "Total tasks: " << total_counter << std::endl;
            std::cout << "Imbalance: " << max_time - min_time << std::endl;
            std::cout << "Imbalance / MAX: " << (max_time - min_time) / max_time * 100 << " %" << std::endl;
            avarege_disb += (max_time - min_time) / max_time * 100;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << tasks_counter << " tasks - rank: " << rank << std::endl;
        std::cout << var_for_calculations << " - current result - rank: " << rank << std::endl;
        update_task();
        tasks_counter = 0;
        iteration_counter++;
    }

    MUTEX_LOCK;
    t1_message_to_send = WORK_DONE;
    MPI_Send(&t1_message_to_send, 1, MPI_INT, rank, REQUEST_TAG, MPI_COMM_WORLD);
    MUTEX_UNLOCK;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) std::cout << "Avg imbalance: " << avarege_disb / WORKING_ITERATIONS << " % " << std::endl;
    return NULL;
}

void* send_task(void* args)
{
    MPI_Status receiver;

    while (iteration_counter < WORKING_ITERATIONS)
    {
        MPI_Recv(&t2_message_to_recv, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &receiver);
        if (receiver.MPI_SOURCE != rank && t2_message_to_recv == NEED_NEW_TASK)
        {
            if (current_own_task >= (rank + 1) * TASKS_IN_LIST - count_of_given_tasks - 1)
            {
                t2_message_to_send = NO_TASK;
                MPI_Send(&t2_message_to_send, 1, MPI_INT, receiver.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
                continue;
            }

            else
            {
                int* task_for_receiver;
                int extra_tasks = 0;
                MUTEX_LOCK;
                count_of_given_tasks += TASKS_TO_GIVE;
                if (count_of_given_tasks + current_own_task >= (rank + 1) * TASKS_IN_LIST - 1)
                {
                    extra_tasks = count_of_given_tasks + current_own_task - (rank + 1) * TASKS_IN_LIST + 1;
                    count_of_given_tasks = (rank + 1) * TASKS_IN_LIST - 1 - current_own_task;
                }
                task_for_receiver = &tasklist[(rank + 1) * TASKS_IN_LIST - count_of_given_tasks - 1];
                MUTEX_UNLOCK;
                t2_message_to_send = SEND_TASK;
                int tasks_to_send = TASKS_TO_GIVE - extra_tasks;
                MPI_Send(&t2_message_to_send, 1, MPI_INT, receiver.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
                MPI_Send(&tasks_to_send, 1, MPI_INT, receiver.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
                MPI_Send(task_for_receiver, tasks_to_send, MPI_INT, receiver.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
                continue;
            }
        }
        else if (t2_message_to_recv == TASK_RECEIVED)
        {
            continue;
        }
        else if (t2_message_to_recv == WORK_DONE)
        {
            break;
        }
    }
    return NULL;
}


int main(int argc, char** argv)
{
    int provided;
    pthread_attr_t attrs;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE)
    {
        std::cout << "Can't init thread" << std::endl;
        MPI_Finalize();
        return -1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    tasklist = new int[TASKS_IN_LIST * proc_num]();

    int check = 0;
    check = pthread_mutex_init(&mutex, NULL);
    if (check != 0)
    {
        std::cout << "Can't init mutex" << std::endl;
        MPI_Finalize();
        return -1;
    }

    check = pthread_attr_init(&attrs);
    if (check != 0)
    {
        std::cout << "Can't init attrs" << std::endl;
        MPI_Finalize();
        return -1;
    }

    check = pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);
    if (check != 0)
    {
        std::cout << "Can't init attrs" << std::endl;
        MPI_Finalize();
        return -1;
    }

    pthread_t threads[2];
    pthread_create(&threads[0], &attrs, send_task, NULL);
    pthread_create(&threads[1], &attrs, resolve_task, NULL);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);
    pthread_attr_destroy(&attrs);
    pthread_mutex_destroy(&mutex);

    MPI_Finalize();
    return 0;
}
