#include <stdlib.h>
#include <mpi.h>
#include <cstring>
#include <assert.h>
#include <time.h>
#include "PP_Labs.cpp"

int main(int argc, char* argv[])
{
    int rank; // ранг процесу
    int size; // загальна кількість процесів

    double tStart;
    double tEnd;
    double tTotal;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N;
    int M;
    bool printMatrix;
    bool verifySolution;
    
    double* matrixAug = NULL;
    double* matrixAugP = NULL;
    double* matrixAugS = NULL;

    if (rank == 0) { // ініціалізація (виконує нульовий процес)
        std::cout << "Size = " << size << "\n";
        InputData(N, M, printMatrix, verifySolution, size);
        std::cout << "\n";
        matrixAug = new double[N * M];
        InitAugmentedMatrix(matrixAug, N);
        matrixAugP = new double[N * M];
        CopyMatrix(matrixAug, matrixAugP, N, M);
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD); // відправляємо/отримуємо введений розмір матриці
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int blockN = N / size; // кількість рядків, що приходиться на ранг
    double* submatrix = new double[blockN * M]; // підматриця для кожного процесу
    MPI_Scatter(matrixAugP, blockN * M, MPI_DOUBLE, submatrix, blockN * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* row = new double[M];

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) { // час початку
        tStart = MPI_Wtime();
    }

    // прямий хід
    double scale;
    for (int i = 0; i < rank * blockN; i ++) { // отримання
        MPI_Bcast(row, M, MPI_DOUBLE, i / blockN, MPI_COMM_WORLD);
        for (int j = 0; j < blockN; j++) { // перетворення рядків
            scale = submatrix[j * M + i];
            for (int k = i + 1; k < M; k++) {
                submatrix[j * M + k] -= scale * row[k];
            }
            submatrix[j * M + i] = 0.0;
        }
    }

    double pivot;
    int i0;
    for (int i = 0; i < blockN; i++) { // відправка
        i0 = rank * blockN + i;
        pivot = submatrix[i * M + i0];
        for (int j = i0 + 1; j < M; j++) { // номралізуємо рядок
            submatrix[i * M + j] /= pivot;
        }
        submatrix[i * M + i0] = 1.0;
        memcpy(row, &submatrix[i * M], M * sizeof(double)); // копіюємо рядок в буфер для відправки
        MPI_Bcast(row, M, MPI_DOUBLE, rank, MPI_COMM_WORLD); // передаємо рядок всім рангам
        for (int j = i + 1; j < blockN; j++) { // перетворення рядків
            scale = submatrix[j * M + i0];
            for (int k = i0 + 1; k < M; k++) {
                submatrix[j * M + k] -= scale * row[k];
            }
            submatrix[j * M + i0] = 0.0;
        }
    }

    for (int i = rank * blockN + blockN; i < N; i++) {
        MPI_Bcast(row, M, MPI_DOUBLE, i / blockN, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD); // синхронізуємо процеси
    
    // зворотній хід
    for (int i = N - 1; i >= rank * blockN + blockN; i--) { // отримання
        MPI_Bcast(row, M, MPI_DOUBLE, i / blockN, MPI_COMM_WORLD);
        for (int j = blockN-1; j >= 0; j--) { // перетворення рядків
            scale = submatrix[j * M + i];
            for (int k = i + 1; k < M; k++) {
                submatrix[j * M + k] -= scale * row[k];
            }
            submatrix[j * M + i] = 0.0;
        }
    }

    for (int i = blockN - 1; i >= 0; i--) { // відправка
        i0 = rank * blockN + i;
        memcpy(row, &submatrix[i * M], M * sizeof(double)); // копіюємо рядок в буфер для відправки
        MPI_Bcast(row, M, MPI_DOUBLE, rank, MPI_COMM_WORLD); // передаємо рядок всім процесам
        for (int j = i - 1; j >= 0; j--) { // перетворення рядків
            scale = submatrix[j * M + i0];
            for (int k = i0 + 1; k < M; k++) {
                submatrix[j * M + k] -= scale * row[k];
            }
            submatrix[j * M + i0] = 0.0;
        }
    }

    for (int i = rank * blockN - 1; i >= 0; i--) {
        MPI_Bcast(row, M, MPI_DOUBLE, i / blockN, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD); // синхронізуємо перед підрахунком часу
    if (rank == 0) { // час завершення
        tEnd = MPI_Wtime();
        tTotal = tEnd - tStart;
    }

    MPI_Gather(submatrix, blockN * M, MPI_DOUBLE, matrixAugP, blockN * M, MPI_DOUBLE, 0, MPI_COMM_WORLD); // збираємо підматриці
    MPI_Finalize(); // завершення

    double* matrixInverseP = NULL;
    double* matrixInverseS = NULL;
    double* matrixOriginal = NULL;

    if (rank == 0) { // вивід результатів
        matrixOriginal = new double[N * N];
        CopyOriginalMatrix(matrixAug, matrixOriginal, N);
        matrixInverseP = new double[N * N];
        CopyInverseMatrix(matrixAugP, matrixInverseP, N);
        
        if (printMatrix) {
            std::cout << "Original:\n";
            PrintMatrix(matrixOriginal, N, N);
            std::cout << "Inverse (parallel):\n";
            PrintMatrix(matrixInverseP, N, N);
        }
        
        if (verifySolution) { // порівняння з послідовним алгоритмом
            matrixAugS = new double[N * M];
            CopyMatrix(matrixAug, matrixAugS, N, M);
            FindInverseSerial(matrixAugS, N);
            matrixInverseS = new double[N * N];
            CopyInverseMatrix(matrixAugS, matrixInverseS, N);
            if (printMatrix) {
                std::cout << "Inverse (parallel):\n";
                PrintMatrix(matrixInverseS, N, N);
            }
            std::cout << "Comparing parallel and serial solutions:\n";
            CompareMatrices(matrixInverseP, matrixInverseS, N, N);
        }

        std::cout << "Time = " << tTotal << "s\n";
    }

    if (rank == 0) { // видалення пам'яті
        delete[] matrixAug;
        delete[] matrixAugP;
        delete[] matrixOriginal;
        delete[] matrixInverseP;
        if (verifySolution) {
            delete[] matrixAugS;
            delete[] matrixInverseS;
        }
    }
    delete[] submatrix;
    delete[] row;
    return 0;
}
