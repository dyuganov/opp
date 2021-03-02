#include <random>
#include <iostream>
#include <time.h>

#include "Matrix.h"

double* nonMpiNonlinearConjugateGradient(double* A, double* b, double* x) {
	double* r = new double[N];
	double* z = new double[N];
	double alpha = 0, beta = 0;
	double* Az = new double[N];
	double* alphaz = new double[N];
	double* betaz = new double[N];
	double* prev_r = new double[N];

	initRandMatrix(A);
	initRandVector(b);
	initRandVector(x);

	double* tmp = new double[N];
	mulMatrixAndVector(A, x, tmp);
	subVector(b, tmp, r);
	delete[] tmp;

	memcpy(z, r, sizeof(double) * N);

	double prevDoProductR = dotProduct(r, r);
	double currDoProductR;
	const double bVectorLength = getVectorLength(b);
	double rVectorLength = getVectorLength(r);

	size_t cnt = 0;
	const double epsilon = 0.00001;

	while (rVectorLength / bVectorLength >= epsilon) {
		mulMatrixAndVector(A, z, Az);
		alpha = dotProduct(r, r) / dotProduct(Az, z);

		mulVectorScalar(z, alpha, alphaz);

		sumVector(x, alphaz, x);
		mulVectorScalar(Az, alpha, Az);

		memcpy(prev_r, r, sizeof(double) * N);

		subVector(r, Az, r);

		currDoProductR = dotProduct(r, r);
		beta = currDoProductR / prevDoProductR;
		prevDoProductR = currDoProductR;

		mulVectorScalar(z, beta, betaz);
		sumVector(r, betaz, z);

		rVectorLength = getVectorLength(r);
		++cnt;
	}

	std::cout << "Iterations: " << cnt << std::endl;

	delete[] alphaz;
	delete[] Az;
	delete[] betaz;
	delete[] prev_r;
	delete[] r;
	delete[] z;

	return x;
}


int main(int argc, char* argv[]) {

	double* A = new double[N * N];
	double* b = new double[N];
	double* x = new double[N];

	clock_t start, end;

	start = clock();
	nonMpiNonlinearConjugateGradient(A, b, x);
	end = clock();

	std::cout << "Time: " << (double)((end - start) / CLOCKS_PER_SEC) << " sec." << std::endl;

	delete[] A;
	delete[] b;
	delete[] x;

	return 0;
}
