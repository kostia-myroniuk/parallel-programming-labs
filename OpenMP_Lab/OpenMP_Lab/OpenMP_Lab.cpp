#include <chrono>
#include <omp.h>
#include "PP_Labs.cpp"

using namespace std::chrono;

void FindInverseOpenMP(double* matrix, int N)
{
    int M = N * 2;
    for (int i = 0; i < N; i++) { 
        #pragma omp parallel for // використовуємо розпаралелювання для внутрішніх циклів
        for (int j = i + 1; j < M; j++) {
            matrix[i * M + j] /= matrix[i * M + i];
        }
        matrix[i * M + i] = 1.0;
        double scale;
        #pragma omp parallel for // використовуємо розпаралелювання для внутрішніх циклів
        for (int j = i + 1; j < N; j++) { 
            for (int k = i + 1; k < M; k++) {
                matrix[j * M + k] -= matrix[i * M + k] * matrix[j * M + i];
            }
            matrix[j * M + i] = 0.0;
        }
    }
    for (int i = N - 1; i >= 0; i--) { 
        double scale;
        #pragma omp parallel for // використовуємо розпаралелювання для внутрішніх циклів
        for (int j = i - 1; j >= 0; j--) { 
            for (int k = i + 1; k < M; k++) {
                matrix[j * M + k] -= matrix[i * M + k] * matrix[j * M + i];
            }
            matrix[j * M + i] = 0.0;
        }
    }
}

int main()
{
    omp_set_num_threads(4);

    int N;
    int M;
    bool printMatrix;
    bool verifySolution;

    InputData(N, M, printMatrix, verifySolution);

    double* matrixAug = new double[N * M];
    InitAugmentedMatrix(matrixAug, N);
    double* matrixAugP = new double[N * M];
    CopyMatrix(matrixAug, matrixAugP, N, M);

    time_point<high_resolution_clock> tStart = high_resolution_clock::now();
    FindInverseOpenMP(matrixAugP, N);
    double tTotal = duration_cast<duration<double, std::ratio<1>>>(high_resolution_clock::now() - tStart).count();

    double* matrixOriginal = new double[N * N];
    CopyOriginalMatrix(matrixAug, matrixOriginal, N);
    double* matrixInverseP = new double[N * N];
    CopyInverseMatrix(matrixAugP, matrixInverseP, N);

    std::cout << "\n";
    if (printMatrix) {
        std::cout << "Original:\n";
        PrintMatrix(matrixOriginal, N, N);
        std::cout << "Inverse (parallel):\n";
        PrintMatrix(matrixInverseP, N, N);
    }
    
    double* matrixAugS = NULL;
    double* matrixInverseS = NULL;
    if (verifySolution) {
        matrixAugS = new double[N * M];
        CopyMatrix(matrixAug, matrixAugS, N, M);
        FindInverseSerial(matrixAugS, N);
        double* matrixInverseS = new double[N * N];
        CopyInverseMatrix(matrixAugS, matrixInverseS, N);
        if (printMatrix) {
            std::cout << "Inverse (serial):\n";
            PrintMatrix(matrixInverseS, N, N);
        }
        std::cout << "Comparing parallel and serial solutions:\n";
        CompareMatrices(matrixInverseP, matrixInverseS, N, N);
    }

    std::cout << "Time = " << tTotal << "s";

    delete[] matrixAug;
    delete[] matrixOriginal;
    delete[] matrixAugP;
    delete[] matrixInverseP;
    if (verifySolution) {
        delete[] matrixAugS;
        delete[] matrixInverseS;
    }
}
