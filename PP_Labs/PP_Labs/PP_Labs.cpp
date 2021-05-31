#include <iostream>
#include <iomanip>

void FindInverseSerial(double* matrix, int N)
{
    int M = N * 2;
    double pivot;
    for (int i = 0; i < N; i++) { // прямий хід
        pivot = matrix[i * M + i]; // за основу беремо елемент на діагоналі матриці А
        for (int j = i + 1; j < M; j++) {
            matrix[i * M + j] /= pivot; // нормалізуємо рядок
        }
        matrix[i * M + i] = 1.0;
        double scale;
        for (int j = i + 1; j < N; j++) { // перетворення рядків
            scale = matrix[j * M + i];
            for (int k = i + 1; k < M; k++) {
                matrix[j * M + k] -= matrix[i * M + k] * scale;
            }
            matrix[j * M + i] = 0.0;
        }
    }
    for (int i = N - 1; i >= 0; i--) { // зворотній хід
        double scale;
        for (int j = i - 1; j >= 0; j--) { // перетворення рядків
            scale = matrix[j * M + i];
            for (int k = i + 1; k < M; k++) {
                matrix[j * M + k] -= matrix[i * M + k] * scale;
            }
            matrix[j * M + i] = 0.0;
        }
    }
}

void InitAugmentedMatrix(double* matrix, int N)
{
    int M = N * 2;
    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            if (j < N) { // випадкове дробове число між -1000 і 1000
                matrix[i * M + j] = (double(rand()) / double(RAND_MAX)) * (1000.0 - -1000.0) + -1000.0;
            }
            else if (i == j - N) {
                matrix[i * M + j] = 1.0;
            }
            else {
                matrix[i * M + j] = 0.0;
            }
        }
    }
}

void CopyMatrix(double* matrixFrom, double* matrixTo, int N, int M)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            matrixTo[i * M + j] = matrixFrom[i * M + j];
        }
    }
}

void CopyOriginalMatrix(double* matrixAugmented, double* matrixTo, int N)
{
    int M = N * 2;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrixTo[i * N + j] = matrixAugmented[i * M + j];
        }
    }
}

void CopyInverseMatrix(double* matrixAugmented, double* matrixTo, int N)
{
    int M = N * 2;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrixTo[i * N + j] = matrixAugmented[i * M + (j + N)];
        }
    }
}

void CompareMatrices(double* matrix1, double* matrix2, int N, int M)
{
    double maxDif = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            if (abs(matrix1[i * M + j] - matrix2[i * M + j]) > maxDif) {
                maxDif = abs(matrix1[i * M + j] - matrix2[i * M + j]);
            }
        }
    }
    if (maxDif < 0.1) {
        std::cout << "Identical\n";
    }
    else {
        std::cout << "Not identical. Max difference = " << maxDif << "\n";
    }
}

void PrintMatrix(double* matrix, int N, int M)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            std::cout << std::setprecision(5) << matrix[i * M + j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void InputData(int& N, int& M, bool& printMatrix, bool& verifySolution, int size = 1)
{
    while (true) {
        std::cout << "Enter N: ";
        std::cin >> N;
        if (N % size == 0) {
            break;
        }
        std::cout << "N must be divisible by size\n";
    }
    M = N * 2;
    int temp;
    std::cout << "Print result (1/0): ";
    std::cin >> temp;
    printMatrix = temp == 1;
    std::cout << "Verify solution (1/0): ";
    std::cin >> temp;
    verifySolution = temp == 1;
}
