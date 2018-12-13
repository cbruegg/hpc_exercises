#include<stdio.h>
#include<stdlib.h>
#include <random>
#include <string>
#include <iostream>
#include <memory>
#include <iomanip>
#include <vector>
#include <sstream>

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
        const auto matrix = generateMatrix(systemSize);
        const auto rhs = generateRhsVector(matrix);

        return 0;
    }

private:

    Matrix generateMatrix(const size_t systemSize) {
        uniform_real_distribution<double> uniform(-1, +1);
        default_random_engine engine;

        auto matrix = vector<vector<double>>(systemSize);
        for (size_t row = 0; row < systemSize; row++) {
            matrix[row] = vector<double>(systemSize);
            for (size_t col = 0; col < systemSize; col++) {
                matrix[row][col] = row == col ? systemSize / 10 : engine();
            }
        }
        return matrix;
    }

    vector<double> generateRhsVector(Matrix matrix) {
        const auto systemSize = matrix.size();
        auto rhs = vector<double>(systemSize);
        for (size_t i = 0; i < systemSize; i++) {
            rhs[i] = matrix[i][i];
        }
        return rhs;
    }

};

int main(int argc, const char *const *const argv) {
    return Main().main(argc, argv);
}