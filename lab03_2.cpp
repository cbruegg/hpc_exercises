#include <vector>
#include <memory>
#include <mpi.h>
#include <utility>
#include <iostream>
#include <random>
#include <chrono>

using namespace std;

int randomStart(int rank) {
    mt19937 engine(static_cast<unsigned long>(chrono::system_clock::now().time_since_epoch().count() - rank * 100));
    uniform_int_distribution<int> uniform(0, 100);
    return uniform(engine);
}

int manualSum(int rank, int totalRanks, int initialValue) {
    const auto phases = log2(totalRanks);
    if (static_cast<int>(phases) != phases) {
        cerr << "Number of machines is not a power of two! Exiting." << endl;
        exit(EXIT_FAILURE);
    }

    auto sum = initialValue;
    for (auto phase = 0; phase < phases; phase++) {
        const auto unitSize = static_cast<int>(pow(2, phase));
        const auto thisUnit = rank / unitSize;
        const auto partner = (thisUnit % 2 == 0) ? rank + unitSize : rank - unitSize;

        int otherSum;
        MPI_Sendrecv(&sum, 1, MPI_INT, partner, 0,
                     &otherSum, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        sum += otherSum;
    }
    return sum;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int totalRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);

    const auto rand = randomStart(rank);
    const auto manualSumResult = manualSum(rank, totalRanks, rand);

    int correctSumResult;
    MPI_Allreduce(&rand, &correctSumResult, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Manual sum computation yielded result " << manualSumResult << endl;
        cout << "Automatic sum computation yielded result " << correctSumResult << endl;

        if (manualSumResult != correctSumResult) {
            cerr << "Result is incorrect!" << endl;
        } else {
            cout << "Result is correct." << endl;
        }
    }

    MPI_Finalize();
    return 0;
}