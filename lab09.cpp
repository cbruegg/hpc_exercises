#include <utility>

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
#include <mpi.h>

#endif

#define DEBUG

using namespace std;

using Matrix = vector<vector<double>>;

int myRank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int ranks() {
    int totalRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);
    return totalRanks;
}

unsigned int rowRemainder(unsigned int systemSize) {
    return systemSize % ranks();
}

unsigned int rowsPerRank(unsigned int systemSize) {
    return systemSize / ranks();
}

class LocalMatrix final {
public:

    const Matrix m;
    const unsigned int rowStart;
    const unsigned int rowEnd;

    LocalMatrix(Matrix m, const unsigned int rowStart, const unsigned int rowEnd) : m(std::move(m)),
                                                                                    rowStart(rowStart),
                                                                                    rowEnd(rowEnd) {}
};

class MatrixOps final {
public:

    static vector<double> times(const LocalMatrix &local, const vector<double> &a) {
        auto localResult = times(local.m, a, local.rowStart, local.rowEnd);
        const auto systemSize = static_cast<unsigned int>(a.size());
        const auto perRank = rowsPerRank(systemSize);
        const auto remainder = rowRemainder(systemSize);

        const auto totalRanks = ranks();

        MPI_Bcast(localResult.data(), perRank + remainder, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (auto bcastRank = 1; bcastRank < totalRanks; bcastRank++) {
            MPI_Bcast(localResult.data() + remainder + bcastRank * perRank,
                      perRank, MPI_DOUBLE, bcastRank, MPI_COMM_WORLD);
        }

        return localResult;
    }

    static vector<double> times(const Matrix &m, const vector<double> &a, unsigned int rowStart, unsigned int rowEnd) {
#ifdef DEBUG
        if (m.size() != a.size() || m[rowStart].size() != a.size()) {
            throw invalid_argument("Sizes are not equal");
        }
#endif

        const auto systemSize = a.size();
        vector<double> c(systemSize);
        for (auto i = rowStart; i < rowEnd; i++) {
            auto sum = 0.0;
            for (auto j = 0u; j < systemSize; j++) {
                sum += m[i][j] * a[j];
            }
            c[i] = sum;
        }

        return c;
    }

    static vector<double> minus(const vector<double> &a, const vector<double> &b) {
#ifdef DEBUG
        if (a.size() != b.size()) {
            throw invalid_argument("Vector sizes are not equal");
        }
#endif

        const auto systemSize = a.size();
        auto c = a;
        for (auto i = 0u; i < systemSize; i++) {
            c[i] -= b[i];
        }
        return c;
    }

    static vector<double> plus(const vector<double> &a, const vector<double> &b) {
#ifdef DEBUG
        if (a.size() != b.size()) {
            throw invalid_argument("Vector sizes are not equal");
        }
#endif

        const auto systemSize = a.size();
        auto c = a;
        for (auto i = 0u; i < systemSize; i++) {
            c[i] += b[i];
        }
        return c;
    }

    static double sqlength(const vector<double> a) {
        auto sum = 0.0;
        for (const auto elem : a) {
            sum += elem * elem;
        }
        return sum;
    }

    static double transposeLeftAndTimes(const vector<double> &a, const vector<double> &b) {
#ifdef DEBUG
        if (a.size() != b.size()) {
            throw invalid_argument("Vector sizes are not equal");
        }
#endif

        const auto systemSize = a.size();
        auto sum = 0.0;
        for (auto i = 0u; i < systemSize; i++) {
            sum += a[i] * b[i];
        }

        return sum;
    }

    static vector<double> times(const double a, const vector<double> &b) {
        auto c = b;
        for (double &i : c) {
            i *= a;
        }
        return c;
    }

    static void println(const Matrix &m) {
        for (auto i = 0u; i < m.size(); i++) {
            for (auto j = 0u; j < m.size(); j++) {
                cout << setw(5) << m[i][j] << " ";
            }
            cout << endl;
        }
    }

    static void println(const LocalMatrix &m) {
        const auto systemSize = m.m[m.rowStart].size();
        for (auto i = 0u; i < systemSize; i++) {
            for (auto j = 0u; j < systemSize; j++) {
                double value;
                if (i < m.rowStart || i >= m.rowEnd) {
                    value = 1 / 0.0;
                } else {
                    value = m.m[i][j];
                }
                cout << setw(5) << value << " ";
            }
            cout << endl;
        }
    }

};

class Main final {
public:
    int main(int argc, char **argv) {
        MPI_Init(&argc, &argv);

        const auto systemSize = stoi(argv[1]);
        const auto sigma = 0.6; // TODO

        if (systemSize <= 0) {
            cerr << "System size needs to be > 0";
            return 1;
        }

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        Matrix matrix;
        vector<double> realSolution;
        vector<double> solution;

        if (rank == 0) {
            matrix = generateMatrix(static_cast<unsigned int>(systemSize), sigma);
            realSolution = generateSolutionVector(static_cast<unsigned int>(systemSize));
            const auto b = generateB(matrix, realSolution);

            solution = solve(&matrix, &b);
        } else {
            solution = solve(nullptr, nullptr);
        }

        if (rank == 0) {
            auto maxErr = 0.0;
            for (auto i = 0u; i < realSolution.size(); i++) {
                maxErr = max(maxErr, abs(solution[i] - realSolution[i]));
            }

            cout << "Max error: " << maxErr << endl;
        }

        MPI_Finalize();
        return 0;
    }

private:

    const double toll = 10e-9;

    LocalMatrix obtainLocalMatrix(const Matrix *const m) {
        const auto rank = myRank();
        const auto totalRanks = ranks();

        unsigned int systemSize;
        if (rank == 0) {
            systemSize = static_cast<int>(m->size());
        }
        MPI_Bcast(&systemSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        const auto perRank = rowsPerRank(systemSize);
        const auto remainder = rowRemainder(systemSize);
        const auto rank0Rows = perRank + remainder;

        vector<double> contigPart;
        if (rank == 0) {
            const auto contig = matrixToContiguousVector(*m);
            contigPart = vector<double>(rank0Rows * systemSize);
            copy(contig.data(), contig.data() + rank0Rows * systemSize, contigPart.data());

            for (auto targetRank = 1; targetRank < totalRanks; targetRank++) {
                MPI_Send(contig.data() + remainder + targetRank * perRank * systemSize, perRank * systemSize,
                         MPI_DOUBLE, targetRank, 0, MPI_COMM_WORLD);
            }
        } else {
            contigPart = vector<double>(perRank * systemSize);
            MPI_Recv(contigPart.data(), perRank * systemSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, nullptr);
        }

        unsigned int rowStart;
        unsigned int rowEnd;
        if (rank == 0) {
            rowStart = 0;
            rowEnd = rank0Rows;
        } else {
            rowStart = remainder + rank * perRank;
            rowEnd = rowStart + perRank;
        }

        Matrix localMatrixPart(systemSize);
        for (auto rowIdx = rowStart; rowIdx < rowEnd; rowIdx++) {
            localMatrixPart[rowIdx] = vector<double>(systemSize);

            copy(contigPart.data() + (rowIdx - rowStart) * systemSize,
                 contigPart.data() + (rowIdx - rowStart + 1) * systemSize,
                 localMatrixPart[rowIdx].data());
        }

        return LocalMatrix(localMatrixPart, rowStart, rowEnd);
    }

    vector<double> matrixToContiguousVector(const Matrix &m) {
        const auto systemSize = m.size();
        const auto totalSize = systemSize * systemSize;
        vector<double> contig(totalSize);
        for (auto i = 0u; i < systemSize; i++) {
            for (auto j = 0u; j < systemSize; j++) {
                contig[i * systemSize + j] = m[i][j];
            }
        }
        return contig;
    }

    vector<double> obtainLocalB(const vector<double> *const b, const unsigned int systemSize) {
        vector<double> localB;

        if (myRank() == 0) {
            localB = *b;
        } else {
            localB = vector<double>(systemSize);
        }

        MPI_Bcast(localB.data(), systemSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        return localB;
    }

    vector<double> solve(const Matrix *const m, const vector<double> *b) {
        auto localMatrix = obtainLocalMatrix(m);
        auto localB = obtainLocalB(b, static_cast<unsigned int>(localMatrix.m.size()));

        vector<double> xk(localB.size(), 0.0);
        auto rk = MatrixOps::minus(localB, MatrixOps::times(localMatrix, xk));
        auto pk = rk;

        const auto bSqLength = MatrixOps::sqlength(localB);
        auto error = numeric_limits<double>::infinity();
        while (error > toll) {
            const auto t = MatrixOps::times(localMatrix, pk);
            const auto alphak = MatrixOps::sqlength(rk) / MatrixOps::transposeLeftAndTimes(pk, t);
            xk = MatrixOps::plus(xk, MatrixOps::times(alphak, pk));
            const auto rk1 = MatrixOps::minus(rk, MatrixOps::times(alphak, t));
            const auto betak = MatrixOps::sqlength(rk1) / MatrixOps::sqlength(rk);
            pk = MatrixOps::plus(rk1, MatrixOps::times(betak, pk));
            rk = rk1;

            error = MatrixOps::sqlength(rk) / bSqLength;
        }

        return xk;
    }

    Matrix generateMatrix(const unsigned int systemSize, const double sigma) {
        const auto s = -0.5;
        const auto d = 1 + sigma;
        Matrix matrix(systemSize);
        for (auto i = 0l; i < systemSize; i++) {
            matrix[i] = vector<double>(systemSize);
            for (auto j = 0l; j < systemSize; j++) {
                if (i == j) {
                    matrix[i][j] = d;
                } else if (abs(i - j) == 1) {
                    matrix[i][j] = s;
                } else {
                    matrix[i][j] = 0;
                }
            }
        }
        matrix[systemSize - 1][0] = s;
        matrix[0][systemSize - 1] = s;

        return matrix;
    }

    vector<double> generateSolutionVector(const unsigned int systemSize) {
        // Have to add epsilon as the distribution generates values from [a, b)
        uniform_real_distribution<double> uniform(0, 1 + numeric_limits<double>::epsilon());
        default_random_engine engine(
                static_cast<unsigned long>(chrono::system_clock::now().time_since_epoch().count())
        );

        auto sol = vector<double>(systemSize);
        for (auto i = 0l; i < systemSize; i++) {
            sol[i] = uniform(engine);
        }

        return sol;
    }

    vector<double> generateB(const Matrix &matrix, const vector<double> &sol) {
        auto b = vector<double>(sol.size());

        for (auto i = 0u; i < b.size(); i++) {
            auto sum = 0.0;
            for (auto j = 0u; j < b.size(); j++) {
                sum += matrix[i][j] * sol[j];
            }
            b[i] = sum;
        }

        return b;
    }

};

int main(int argc, char **argv) {
    return Main().main(argc, argv);
}