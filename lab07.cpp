#include<stdio.h>
#include<stdlib.h>
#include <random>
#include <string>
#include <iostream>
#include <memory>
#include <iomanip>
#include <vector>
#include <sstream>
#include <optional>

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

        printSystem(matrix, rhs);
        cout << endl;
        gaussianElimination(matrix, rhs);
        printSystem(matrix, rhs);
        cout << endl;

        const auto solution = solveWithBackSubstitution(matrix, rhs);
        for (auto value : solution) {
            cout << setw(15) << value << " ";
        }
        cout << endl;

        return 0;
    }

private:

    Matrix generateMatrix(const size_t systemSize) {
        uniform_real_distribution<double> uniform(0, 1);
        default_random_engine engine;

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
            for (size_t row = 0; row < col; row++) {
                sol[row] -= matrix[row][col] * sol[col];
            }

            if (col-- == 0) break;
        }

        return sol;
    }

    void printSystem(Matrix matrix, vector<double> rhs) {
        for (size_t rowIdx = 0; rowIdx < matrix.size(); rowIdx++) {
            const auto row = matrix[rowIdx];
            for (const auto cell : row) {
                cout << setw(15) << cell << " ";
            }
            cout << setw(15) << rhs[rowIdx] << " ";
            cout << endl;
        }
    }

};

int main(int argc, const char *const *const argv) {
    return Main().main(argc, argv);
}