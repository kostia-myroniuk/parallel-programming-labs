#include <chrono>
#include "PP_Labs.cpp"

void InitIdentityMatrix(double* matrixI, int N)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrixI[i * N + j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

void MultiplySquareMatrices(double* result, double* A, double* B, int N)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i * N + j] = 0.0;
            for (int k = 0; k < N; k++) {
                result[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}

void VerifyMultiplication(double* matrixOriginal, double* matrixResult, int N, bool& printMatrix)
{
    double* matrixI = new double[N * N];
    InitIdentityMatrix(matrixI, N);
    double* matrixMult = new double[N * N];

    MultiplySquareMatrices(matrixMult, matrixOriginal, matrixResult, N);
    std::cout << "Multiplication (A*B)\n";
    std::cout << "Comparing with identity matrix:\n";
    CompareMatrices(matrixMult, matrixI, N, N);
    if (printMatrix) {
        PrintMatrix(matrixMult, N, N);
    }

    MultiplySquareMatrices(matrixMult, matrixResult, matrixOriginal, N);
    std::cout << "Multiplication (B*A)\n";
    std::cout << "Comparing with identity matrix:\n";
    CompareMatrices(matrixMult, matrixI, N, N);
    if (printMatrix) {
        PrintMatrix(matrixMult, N, N);
    }

    delete[] matrixI;
    delete[] matrixMult;
}

int main()
{
    using namespace std::chrono;

    int N;
    int M;
    bool printMatrix;
    bool verifySolution;
    InputData(N, M, printMatrix, verifySolution);

    double* matrixAug = new double[N * M];
    InitAugmentedMatrix(matrixAug, N);

    double* matrixOriginal = new double[N * N];
    CopyOriginalMatrix(matrixAug, matrixOriginal, N);

    time_point<high_resolution_clock> tStart = high_resolution_clock::now();
    FindInverseSerial(matrixAug, N);
    double tTotal = duration_cast<duration<double, std::ratio<1>>>(high_resolution_clock::now() - tStart).count();
    
    double* matrixInverse = new double[N * N];
    CopyInverseMatrix(matrixAug, matrixInverse, N);

    std::cout << "\n";
    if (printMatrix) {
        std::cout << "Original:\n";
        PrintMatrix(matrixOriginal, N, N);
        std::cout << "Inverse:\n";
        PrintMatrix(matrixInverse, N, N);
    }
    
    if (verifySolution) {
        VerifyMultiplication(matrixOriginal, matrixInverse, N, printMatrix);
    }

    std::cout << "Time = " << tTotal << "s";

    delete[] matrixAug;
    delete[] matrixOriginal;
    delete[] matrixInverse;
}
