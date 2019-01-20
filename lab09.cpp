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
#include <mpi.h>

#ifdef _OPENMP

#include <omp.h>

#endif

#define DEBUG false

using namespace std;

using Matrix = vector<vector<double>>;

unsigned int myRank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return static_cast<unsigned int>(rank);
}

unsigned int ranks() {
    int totalRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);
    return static_cast<unsigned int>(totalRanks);
}

unsigned int rowRemainder(unsigned int systemSize) {
    return systemSize % ranks();
}

unsigned int rowsPerRank(unsigned int systemSize) {
    return systemSize / ranks();
}

unsigned int rowStart(unsigned int systemSize, int rank = myRank()) {
    if (rank == 0) {
        return 0;
    } else {
        const auto perRank = rowsPerRank(systemSize);
        const auto remainder = rowRemainder(systemSize);
        return remainder + rank * perRank;
    }
}

unsigned int rowEnd(unsigned int systemSize, int rank = myRank()) {
    const auto perRank = rowsPerRank(systemSize);
    const auto remainder = rowRemainder(systemSize);

    if (rank == 0) {
        return perRank + remainder;
    } else {
        const auto rowStart = remainder + rank * perRank;
        return rowStart + perRank;
    }
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
        return times(local.m, a, local.rowStart, local.rowEnd);
    }

    static vector<double> times(const Matrix &m, const vector<double> &a, unsigned int rowStart, unsigned int rowEnd) {
#ifdef DEBUG
        if (m.size() != a.size() || m[rowStart].size() != a.size()) {
            throw invalid_argument("Sizes are not equal");
        }
#endif

        const auto systemSize = a.size();
        vector<double> c(systemSize);
#pragma omp parallel for
        for (auto i = rowStart; i < rowEnd; i++) {
            auto sum = 0.0;
            for (auto j = 0u; j < systemSize; j++) {
                sum += m[i][j] * a[j];
            }
            c[i] = sum;
        }

        return c;
    }

    static vector<double> localMinus(const vector<double> &a, const vector<double> &b) {
#if DEBUG
        if (a.size() != b.size()) {
            throw invalid_argument("Vector sizes are not equal");
        }
#endif

        const auto systemSize = static_cast<unsigned int>(a.size());
        const auto start = rowStart(systemSize);
        const auto end = rowEnd(systemSize);

        auto c = a;
#pragma omp parallel for
        for (auto i = start; i < end; i++) {
            c[i] -= b[i];
        }
        return c;
    }

    static vector<double> localPlus(const vector<double> &a, const vector<double> &b) {
#if DEBUG
        if (a.size() != b.size()) {
            throw invalid_argument("Vector sizes are not equal");
        }
#endif

        const auto systemSize = static_cast<unsigned int>(a.size());
        const auto start = rowStart(systemSize);
        const auto end = rowEnd(systemSize);

        auto c = a;
#pragma omp parallel for
        for (auto i = start; i < end; i++) {
            c[i] += b[i];
        }
        return c;
    }

    static double sqlength(const vector<double> &a) {
        return transposeLeftAndTimes(a, a);
    }

    static double transposeLeftAndTimes(const vector<double> &a, const vector<double> &b) {
#if DEBUG
        if (a.size() != b.size()) {
            throw invalid_argument("Vector sizes are not equal");
        }
#endif

        const auto systemSize = static_cast<unsigned int>(a.size());
        const auto start = rowStart(systemSize);
        const auto end = rowEnd(systemSize);

        auto localSum = 0.0;
#pragma omp parallel for reduction(+: localSum)
        for (auto i = start; i < end; i++) {
            localSum += a[i] * b[i];
        }

        auto globalSum = localSum;
        MPI_Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return globalSum;
    }

    static vector<double> localTimes(const double a, const vector<double> &b) {
        const auto systemSize = static_cast<unsigned int>(b.size());
        const auto start = rowStart(systemSize);
        const auto end = rowEnd(systemSize);

        auto c = b;
#pragma omp parallel for
        for (auto i = start; i < end; i++) {
            c[i] *= a;
        }
        return c;
    }

};

class Main final {
public:
    int main(int argc, char **argv) {
        MPI_Init(&argc, &argv);

        const auto systemSize = stoi(argv[1]);
        const auto sigma = 0.6;

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
                MPI_Send(contig.data() + remainder * systemSize + targetRank * perRank * systemSize,
                         perRank * systemSize,
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

    vector<double> obtainLocalB(const vector<double> *const b, const unsigned int systemSize,
                                const vector<int> &vecRowRecvCounts, const vector<int> &vecRecvDspls) {
        const double* bData;

        if (myRank() == 0) {
            bData = b->data();
        } else {
            bData = nullptr;
        }

        const auto start = rowStart(systemSize);
        const auto end = rowEnd(systemSize);
        auto localB = vector<double>(systemSize);
        MPI_Scatterv(bData, vecRowRecvCounts.data(), vecRecvDspls.data(), MPI_DOUBLE, localB.data() + start,
                     end - start, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        return localB;
    }

    vector<double> solve(const Matrix *const m, const vector<double> *b) {
        const auto localMatrix = obtainLocalMatrix(m);
        const auto systemSize = static_cast<unsigned int>(localMatrix.m.size());

        auto vecRecvBuf = vector<double>(systemSize);
        auto vecRowRecvCounts = vector<int>(ranks());
        auto vecRecvDispls = vector<int>(ranks());
        for (int i = 0, rowSum = 0; i < vecRowRecvCounts.size(); i++) {
            auto rowCountOfRankI = rowEnd(systemSize, i) - rowStart(systemSize, i);

            vecRowRecvCounts[i] = rowCountOfRankI;
            vecRecvDispls[i] = rowSum;

            rowSum += rowCountOfRankI;
        }

        const auto localB = obtainLocalB(b, static_cast<unsigned int>(localMatrix.m.size()),
                                         vecRowRecvCounts, vecRecvDispls);

        const auto start = rowStart(systemSize);
        const auto end = rowEnd(systemSize);

        vector<double> xk(localB.size(), 0.0);
        auto rk = MatrixOps::localMinus(localB, MatrixOps::times(localMatrix, xk));
        auto pk = rk;

        const auto bSqLength = MatrixOps::sqlength(localB);
        auto error = numeric_limits<double>::infinity();
        while (error > toll) {
            // Redistribute pk
            MPI_Allgatherv(pk.data() + start, end - start, MPI_DOUBLE, vecRecvBuf.data(), vecRowRecvCounts.data(),
                           vecRecvDispls.data(),
                           MPI_DOUBLE,
                           MPI_COMM_WORLD);
            pk = vecRecvBuf;

            // Parallelized. Requires full pk
            const auto t = MatrixOps::times(localMatrix, pk);

            // Parallelized. Only requires local rows
            const auto alphak = MatrixOps::sqlength(rk) / MatrixOps::transposeLeftAndTimes(pk, t);

            // Parallelized. Only requires local rows
            xk = MatrixOps::localPlus(xk, MatrixOps::localTimes(alphak, pk));

            // Parallelized. Only requires local rows
            const auto rk1 = MatrixOps::localMinus(rk, MatrixOps::localTimes(alphak, t));

            // Parallelized. Only requires local rows
            const auto betak = MatrixOps::sqlength(rk1) / MatrixOps::sqlength(rk);

            // Parallelized. Only requires local rows
            pk = MatrixOps::localPlus(rk1, MatrixOps::localTimes(betak, pk));

            // Parallelized. Only requires local rows
            rk = rk1;

            // Parallelized. Only requires local rows
            error = MatrixOps::sqlength(rk) / bSqLength;
        }

        // Redistribute xk
        MPI_Allgatherv(xk.data() + start, end - start, MPI_DOUBLE, vecRecvBuf.data(), vecRowRecvCounts.data(),
                       vecRecvDispls.data(),
                       MPI_DOUBLE,
                       MPI_COMM_WORLD);
        xk = vecRecvBuf;

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