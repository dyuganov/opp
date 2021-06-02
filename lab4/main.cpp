#include <iostream>
#include <cmath>

#ifdef __unix__
#include <mpi.h>
#elif defined(_WIN32) || defined(WIN32)
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#endif

const double eps = 0.00000001;
const int a = 100000;

const double Dx = 2;
const double Dy = 2;
const double Dz = 2;

const int Nx = 320;
const int Ny = 320;
const int Nz = 320;

const double hx = (Dx / (Nx - 1));
const double hy = (Dy / (Ny - 1));
const double hz = (Dz / (Nz - 1));

const double multiplier = 1 / (2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz) + a);

double phi(double x, double y, double z){
    double phi = 0;
    phi = x * x + y * y + z * z;
    return phi;
}

double rho(double x, double y, double z){
    double rho = 0;
    rho = 6 - a * phi(x, y, z);
    return rho;
}

void calculateCenter(double layerHeight, double* prevPhi, double* Phi, int rank, bool* flag){
    double xComp, yComp, zComp;

    for (int k = 1; k < layerHeight - 1; k++){
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < Nx-1; j++){
                xComp = (prevPhi[Nx * Ny * k + Nx * i + (j - 1)] + prevPhi[Nx * Ny * k + Nx * i + (j + 1)]) / (hx * hx);
                yComp = (prevPhi[Nx * Ny * k + Nx * (i - 1) + j] + prevPhi[Nx * Ny * k + Nx * (i + 1) + j]) / (hy * hy);
                zComp = (prevPhi[Nx * Ny * (k - 1) + Nx * i + j] + prevPhi[Nx * Ny * (k + 1) + Nx * i + j]) / (hz * hz);

                Phi[Nx * Ny * k + Nx * i + j] = multiplier * (xComp + yComp + zComp - rho(-Dx / 2 + j * hx, -Dy / 2 + i * hy, -Dz / 2 + (k + layerHeight * rank) * hz));

                if (fabs(Phi[Nx * Ny * k + Nx * i + j] - prevPhi[Nx * Ny * k + Nx * i + j]) > eps){
                    (*flag) = true;
                }
            }
        }
    }
}

void calculateEdges(int layerHeight, double* prevPhi, double* Phi, int rank, bool* flag, const double* downLayer, const double* upLayer, int procNum){
    double xComp, yComp, zComp;

    for (int i = 1; i < Ny - 1; i++){
        for (int j = 1; j < Nx - 1; j++){
            if (rank != 0){
                xComp = (prevPhi[Nx * i + (j - 1)] + prevPhi[Nx * i + (j + 1)]) / (hx * hx);
                yComp = (prevPhi[Nx * (i - 1) + j] + prevPhi[Nx * (i + 1) + j]) / (hy * hy);
                zComp = (downLayer[Nx * i + j] + prevPhi[Nx * Ny + Nx * i + j]) / (hz * hz);

                Phi[Nx * i + j] = multiplier * (xComp + yComp + zComp - rho(-Dx / 2 + j * hx, -Dy / 2 + i * hy, -Dz / 2 + (layerHeight * rank) * hz));

                if (fabs(Phi[Nx * i + j] - prevPhi[Nx * i + j]) > eps){
                    (*flag) = true;
                }
            }

            if (rank != procNum - 1){
                xComp = (prevPhi[Nx * Ny * (layerHeight - 1) + Nx * i + (j - 1)] + prevPhi[Nx * Ny * (layerHeight - 1) + Nx * i + (j + 1)]) / (hx * hx);
                yComp = (prevPhi[Nx * Ny * (layerHeight - 1) + Nx * (i - 1) + j] + prevPhi[Nx * Ny * (layerHeight - 1) + Nx * (i + 1) + j]) / (hy * hy);
                zComp = (prevPhi[Nx * Ny * (layerHeight - 2) + Nx * i + j] + upLayer[Nx * i + j]) / (hz * hz);

                Phi[Nx * Ny * (layerHeight - 1) + Nx * i + j] = multiplier * (xComp + yComp + zComp - rho(-Dx / 2 + j * hx, -Dy / 2 + i * hy, -Dz / 2 + ((layerHeight - 1) + layerHeight * rank) * hz));

                if (fabs(Phi[Nx * i + j] - prevPhi[Nx * i + j]) > eps){
                    (*flag) = true;
                }
            }
        }
    }
}

void calcMaxDifference(int rank, int layerHeight, double* Phi ){
    double max = 0;
    double diff = 0;

    for (int k = 0; k < layerHeight; k++){
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < Nx; j++){
                diff = fabs(Phi[k * Nx * Ny + i * Nx + j] - phi(-Dx / 2 + j * hx, -Dy / 2 + i * hy, -Dz / 2 + (k + layerHeight * rank) * hz));
                if (diff > max){
                    max = diff;
                }
            }
        }
    }

    double tmp = 0;
    MPI_Allreduce(&max, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0){
        max = tmp;
        std::cout << "Max difference: " << max << std::endl;
    }
}

int main(int argc, char** argv){
    int rank, procNum;
    double timeStart, timeFinish;
    bool deltaLargerEps = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);

    int layerHeight = Nz / procNum;
    double* Phi = new double[Nx * Ny * layerHeight];
    double* prevPhi = new double[Nx * Ny * layerHeight];
    double* downLayer = new double[Nx * Ny];
    double* upLayer = new double[Nx * Ny];

    if (rank == 0) timeStart = MPI_Wtime();

    for (int k = 0; k < layerHeight; k++){
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < Nx; j++){
                if (i == 0 || j == 0 || i == Ny - 1 || j == Nx - 1){
                    Phi[Nx * Ny * k + Nx * i + j] = phi(-Dx / 2 + j * hx, -Dy / 2 + i * hy, -Dz/2 + (k + layerHeight * rank) * hz);
                    prevPhi[Nx * Ny * k + Nx * i + j] = phi(-Dx / 2 + j * hx, -Dy / 2 + i * hy, -Dz / 2 + (k + layerHeight * rank) * hz);
                }
                else{
                    Phi[Nx * Ny * k + Nx * i + j] = 0;
                    prevPhi[Nx * Ny * k + Nx * i + j] = 0;
                }
            }
        }
    }
    if (rank == 0){
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < Nx; j++){
                Phi[0 + Nx * i + j] = phi(-Dx / 2 + j * hx, -Dy / 2 + i * hy, -Dz / 2);
                prevPhi[0 + Nx * i + j] = phi(-Dx / 2 + j * hx, -Dy / 2 + i * hy, -Dz / 2);
            }
        }
    }

    if (rank == procNum - 1){
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < Nx; j++){
                Phi[Nx * Ny * (layerHeight - 1) + Nx * i + j] = phi(-Dx / 2 + j * hx, -Dy / 2 + i * hy, Dz / 2);
                prevPhi[Nx * Ny * (layerHeight - 1) + Nx * i + j] = phi(-Dx / 2 + j * hx, -Dy / 2 + i * hy, Dz / 2);
            }
        }
    }

    double* tmp;
    bool lor_res;
    double counter = 0;
    MPI_Request requests[4];

    while (deltaLargerEps){
        deltaLargerEps = false;

        tmp = prevPhi;
        prevPhi = Phi;
        Phi = tmp;

        if (rank != 0){
            MPI_Isend(&prevPhi[0], Nx * Ny, MPI_DOUBLE, rank - 1, 10, MPI_COMM_WORLD, &requests[0]); //MPI_Request
            MPI_Irecv(downLayer, Nx * Ny, MPI_DOUBLE, rank - 1, 20, MPI_COMM_WORLD, &requests[1]);
        }

        if (rank != procNum - 1){
            MPI_Isend(&prevPhi[(layerHeight - 1) * Nx * Ny], Nx * Ny, MPI_DOUBLE, rank + 1, 20, MPI_COMM_WORLD, &requests[2]); //MPI_Request
            MPI_Irecv(upLayer, Nx * Ny, MPI_DOUBLE, rank + 1, 10, MPI_COMM_WORLD, &requests[3]);
        }

        calculateCenter(layerHeight, prevPhi, Phi, rank, &deltaLargerEps);

        if (rank != 0){
            MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
            MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
        }

        if (rank != procNum - 1){
            MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
            MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
        }

        calculateEdges(layerHeight, prevPhi, Phi, rank, &deltaLargerEps, downLayer, upLayer, procNum);

        MPI_Allreduce(&deltaLargerEps, &lor_res, 1, MPI_CHAR, MPI_LOR, MPI_COMM_WORLD);
        deltaLargerEps = lor_res;

        counter++;
    }

    if (rank == 0) timeFinish = MPI_Wtime();

    double tmpCounter = 0;
    MPI_Allreduce(&counter, &tmpCounter, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    counter = tmpCounter;

    calcMaxDifference(rank, layerHeight, Phi);

    if (rank == 0){
        std::cout << "Number of iterations:" << counter << std::endl;
        std::cout << "Time:" << (timeFinish - timeStart) << std::endl;
    }

    MPI_Finalize();
    return 0;
}