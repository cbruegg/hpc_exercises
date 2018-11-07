#include <vector>
#include <memory>
#include <random>
#include <limits>
#include <mpi.h>
#include <utility>
#include <iostream>

using namespace std;

uint64_t tossesInCircle(uint64_t totalTosses) {
    uniform_real_distribution<double> uniform(-1, +1);
    default_random_engine engine;

    uint64_t result = 0;
    for (uint64_t toss = 0; toss < totalTosses; toss++) {
        const auto x = uniform(engine);
        const auto y = uniform(engine);
        const auto distanceSquared = x * x + y * y;
        if (distanceSquared <= 1) {
            result++;
        }
    }

    return result;
}

uint64_t stou64(const string &s) {
    const unsigned long long result = stoull(s);
    if (result > numeric_limits<uint64_t>::max()) {
        throw invalid_argument("Number is too large!");
    }
    return result;
}

uint64_t readTossesFromCinOrExit() {
    cout << "How often should I toss?" << endl;

    string line;
    getline(cin, line);

    try {
        return stou64(line);
    }
    catch (invalid_argument &e) {
        cerr << "Invalid input!" << endl;
        exit(EXIT_FAILURE);
    };
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int totalRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);

    uint64_t tosses;
    if (rank == 0) {
        tosses = readTossesFromCinOrExit();
    }
    MPI_Bcast(&tosses, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    const auto tossesPerRank = tosses / totalRanks;
    const auto tossesInCircleResult = tossesInCircle(tossesPerRank);
    uint64_t globalTossesInCircleResult;
    MPI_Reduce(&tossesInCircleResult, &globalTossesInCircleResult, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        const auto piEstimate = (4.0 * globalTossesInCircleResult) / tosses;
        cout << globalTossesInCircleResult << " tosses landed inside the circle. π ≈ " << piEstimate << endl;
    }

    MPI_Finalize();
    return 0;
}
