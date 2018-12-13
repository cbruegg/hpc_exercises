// Scheduler: static. Chunk size: 3000. Threads: 1.
//   Max error in solution = 2.40918e-14
//   Time for Gaussian elim = 12 seconds
//   Time for back sub = 0 seconds
//   Total time for solve = 12 seconds
// BASELINE

// Scheduler: static. Chunk size: 75. Threads: 40.
//   Max error in solution = 2.95319e-14
//   Time for Gaussian elim = 0 seconds
//   Time for back sub = 0 seconds
//   Total time for solve = 0 seconds
// SPEEDUP > 12

// Scheduler: static. Chunk size: 25. Threads: 40.
//   Max error in solution = 2.39808e-14
//   Time for Gaussian elim = 0 seconds
//   Time for back sub = 0 seconds
//   Total time for solve = 0 seconds
// SPEEDUP > 12

// Scheduler: static. Chunk size: 1. Threads: 40.
//   Max error in solution = 2.81997e-14
//   Time for Gaussian elim = 1 seconds
//   Time for back sub = 0 seconds
//   Total time for solve = 1 seconds
// SPEEDUP ≈ 12

// Scheduler: dynamic. Chunk size: 1. Threads: 40.
//   Max error in solution = 2.90878e-14
//   Time for Gaussian elim = 1 seconds
//   Time for back sub = 0 seconds
//   Total time for solve = 1 seconds
// SPEEDUP ≈ 12

// Scheduler: dynamic. Chunk size: 25. Threads: 40.
//   Max error in solution = 2.94209e-14
//   Time for Gaussian elim = 1 seconds
//   Time for back sub = 0 seconds
//   Total time for solve = 1 seconds
// SPEEDUP ≈ 12

// Scheduler: guided. Chunk size: 1. Threads: 40.
//   Max error in solution = 2.93099e-14
//   Time for Gaussian elim = 1 seconds
//   Time for back sub = 0 seconds
//   Total time for solve = 1 seconds
// SPEEDUP ≈ 12

// Scheduler: guided. Chunk size: 25. Threads: 40.
//   Max error in solution = 3.24185e-14
//   Time for Gaussian elim = 1 seconds
//   Time for back sub = 0 seconds
//   Total time for solve = 1 seconds
// SPEEDUP ≈ 12

// Result: Static scheduling with equally-sized blocks seems to be the most
//         efficient method in this case. This suggests that the workload
//         per block is balanced well. Taking a look at the bodies of the
//         parallelized for-loops confirms this: the amount of work
//         per iteration is mostly dependent on the system size for the
//         gaussian elimination and even constant for the back substitution.

#include<stdio.h>
#include<stdlib.h>
#include <random>
#include <string>
#include <iostream>
#include <memory>
#include <iomanip>
#include <vector>
#include <sstream>
#include <chrono>

#ifdef _OPENMP

#include <omp.h>

#endif

using namespace std;

using Matrix = vector<vector<double>>;

class Main final {
public:
    int main(const int argc, const char *const *const argv) {
        if (argc < 2) {
            cerr << "No system size specified!" << endl;
            exit(1);
        }

        const auto systemSize = static_cast<const size_t>(stoi(argv[1]));
        auto matrix = generateMatrix(systemSize);
        auto rhs = generateRhsVector(matrix);

        const auto elimBegin = chrono::steady_clock::now();
        gaussianElimination(matrix, rhs);
        const auto elimEnd = chrono::steady_clock::now();

        const auto solveBegin = chrono::steady_clock::now();
        const auto solution = solveWithBackSubstitution(matrix, rhs);
        const auto solveEnd = chrono::steady_clock::now();

        auto maxErr = 0.0;
        for (auto value : solution) {
            maxErr = max(maxErr, abs(1.0 - value));
        }
        cout << "Max error in solution = " << maxErr << endl;
        cout << "Time for Gaussian elim = "
             << chrono::duration_cast<chrono::seconds>(elimEnd - elimBegin).count()
             << " seconds"
             << endl;
        cout << "Time for back sub = "
             << chrono::duration_cast<chrono::seconds>(solveEnd - solveBegin).count()
             << " seconds"
             << endl;
        cout << "Total time for solve = "
             << chrono::duration_cast<chrono::seconds>(elimEnd - elimBegin + solveEnd - solveBegin).count()
             << " seconds"
             << endl;


        return 0;
    }

private:

    Matrix generateMatrix(const size_t systemSize) {
        uniform_real_distribution<double> uniform(0, 1);
        default_random_engine engine(
                static_cast<unsigned long>(chrono::system_clock::now().time_since_epoch().count())
        );

        auto matrix = vector<vector<double>>(systemSize);
        for (size_t row = 0; row < systemSize; row++) {
            matrix[row] = vector<double>(systemSize);
            for (size_t col = 0; col < systemSize; col++) {
                matrix[row][col] = row == col ? systemSize / 10.0 : uniform(engine);
            }
        }
        return matrix;
    }

    vector<double> generateRhsVector(Matrix matrix) {
        const auto systemSize = matrix.size();
        auto rhs = vector<double>(systemSize);
        for (size_t row = 0; row < systemSize; row++) {
            rhs[row] = 0;

            for (size_t col = 0; col < systemSize; col++) {
                rhs[row] += matrix[row][col];
            }
        }
        return rhs;
    }

    void gaussianElimination(Matrix &matrix, vector<double> &rhs) {
        const auto systemSize = matrix.size();

        for (size_t col = 0; col < systemSize; col++) {
#pragma omp parallel for schedule(runtime)
            for (auto row = col + 1; row < systemSize; row++) {
                const auto eliminationFactor = matrix[row][col] / matrix[col][col];
                for (size_t i = 0; i < systemSize; i++) {
                    matrix[row][i] -= eliminationFactor * matrix[col][i];
                }
                rhs[row] -= eliminationFactor * rhs[col];
            }
        }
    }

    vector<double> solveWithBackSubstitution(const Matrix &matrix, const vector<double> &rhs) {
        const auto systemSize = matrix.size();
        auto sol = rhs;

        for (auto col = systemSize - 1;;) {
            sol[col] /= matrix[col][col];

#pragma omp parallel for schedule(runtime)
            for (size_t row = 0; row < col; row++) {
                sol[row] -= matrix[row][col] * sol[col];
            }

            if (col-- == 0) break;
        }

        return sol;
    }

};

int main(int argc, const char *const *const argv) {
    return Main().main(argc, argv);
}