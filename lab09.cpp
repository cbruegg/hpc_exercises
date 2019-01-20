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

#define DEBUG

using namespace std;

using Matrix = vector<vector<double>>;

class MatrixOps final {
public:
    static vector<double> times(const Matrix &m, const vector<double> &a) {
#ifdef DEBUG
        if (m.size() != a.size() || m[0].size() != a.size()) {
            throw invalid_argument("Sizes are not equal");
        }
#endif

        const auto systemSize = a.size();
        vector<double> c(systemSize);
        for (auto i = 0u; i < systemSize; i++) {
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

};

class Main final {
public:
    int main(const int argc, const char *const *const argv) {

        const auto systemSize = stoi(argv[1]);
        const auto sigma = 0.6; // TODO

        if (systemSize <= 0) {
            cerr << "System size needs to be > 0";
            return 1;
        }

        const auto matrix = generateMatrix(static_cast<unsigned int>(systemSize), sigma);
        const auto realSolution = generateSolutionVector(static_cast<unsigned int>(systemSize));
        const auto b = generateB(matrix, realSolution);

        const auto solution = solve(matrix, b);

        auto maxErr = 0.0;
        for (auto i = 0u; i < realSolution.size(); i++) {
            maxErr = max(maxErr, abs(solution[i] - realSolution[i]));
        }

        cout << "Max error: " << maxErr << endl;

        return 0;
    }

private:

    const double toll = 10e-9;

    vector<double> solve(const Matrix &m, const vector<double> &b) {
        vector<double> xk(b.size(), 0.0);
        auto rk = MatrixOps::minus(b, MatrixOps::times(m, xk));
        auto pk = rk;

        const auto bSqLength = MatrixOps::sqlength(b);
        auto error = numeric_limits<double>::infinity();
        while (error > toll) {
            const auto t = MatrixOps::times(m, pk);
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

int main(int argc, const char *const *const argv) {
    return Main().main(argc, argv);
}