#include <vector>
#include <memory>
#include <mpi.h>
#include <utility>
#include <limits>
#include <iostream>

using namespace std;

class HelloWorld {

public:
    static void sayHello() {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        cout << "Hello! I'm rank " << rank << " out of " << size << '\n';
    }
};

class PingPong {
private:
    static const int BENCHMARK_REPETITIONS = 50;

    static void sendAndReceive(const bool print = false) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        switch (rank) {
            case 0: {
                auto number = 1;
                MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (print) { cout << "Pong\n"; }
                break;
            }
            case 1: {
                int number;
                MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (print) { cout << "Ping\n"; }
                MPI_Send(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            }
            default:
                break;
        }
    }

public:
    static void warmup() {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (rank == 0) { cout << "\nWarming up latency benchmark...\n"; }
        sendAndReceive(true);
        if (rank == 0) { cout << "Warmup complete.\n"; }
    }

    static void benchmark() {
        const auto startTimeS = MPI_Wtime();
        for (auto i = 0; i < BENCHMARK_REPETITIONS; i++) {
            sendAndReceive(false);
        }
        const auto endTimeS = MPI_Wtime();
        const auto totalTimeMicroS = (endTimeS - startTimeS) * 1e6;
        const auto transmissionTimeMicroS = totalTimeMicroS / BENCHMARK_REPETITIONS;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) { cout << "Time per transmission (back and forth): " << transmissionTimeMicroS << " µs\n\n"; }
    }
};

class BlobBenchmark {
private:
    static const int BENCHMARK_REPETITIONS = 50;

    static void sendAndReceive(shared_ptr<vector<int32_t>> arr, const bool print = false) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int size;
        if (arr->size() > INT32_MAX) {
            throw invalid_argument("Array is too large!");
        } else {
            size = static_cast<int>(arr->size());
        }

        switch (rank) {
            case 0: {
                MPI_Send(arr->data(), size, MPI_INT32_T, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(arr->data(), size, MPI_INT32_T, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (print) { cout << "Pong\n"; }
                break;
            }
            case 1: {
                MPI_Recv(arr->data(), size, MPI_INT32_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (print) { cout << "Ping\n"; }
                MPI_Send(arr->data(), size, MPI_INT32_T, 0, 0, MPI_COMM_WORLD);
                break;
            }
            default:
                break;
        }
    }

public:
    static void warmup(shared_ptr<vector<int32_t>> arr) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (rank == 0) { cout << "Warming up blob benchmark...\n"; }
        sendAndReceive(move(arr), true);
        if (rank == 0) { cout << "Warmup complete.\n"; }
    }

    static void benchmark(shared_ptr<vector<int32_t>> arr) {
        const auto startTimeS = MPI_Wtime();
        for (auto i = 0; i < BENCHMARK_REPETITIONS; i++) {
            sendAndReceive(arr, false);
        }
        const auto endTimeS = MPI_Wtime();
        const auto totalTimeS = endTimeS - startTimeS;
        const auto transmissionTimeS = totalTimeS / (BENCHMARK_REPETITIONS * 2);
        const auto gbPerS = transmissionTimeS == 0.0 ?
                            std::numeric_limits<double>::infinity() :
                            (arr->size() * 8.0 / 1000000000) / transmissionTimeS;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) { cout << "One-way transmission speed: " << gbPerS << " Gb/s\n"; }
    }
};

const auto MIB = 1048576UL /* 1 MiB */;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    HelloWorld::sayHello();
    MPI_Barrier(MPI_COMM_WORLD);

    PingPong::warmup();
    PingPong::benchmark();

    for (const auto benchmarkSize : {1 * MIB, 2 * MIB, 4 * MIB, 16 * MIB,
                                     32 * MIB}) {
        auto arr = make_shared<vector<int32_t >>(vector<int32_t>());
        arr->resize(benchmarkSize / 8);
        std::fill(arr->begin(), arr->end(), 0);
        BlobBenchmark::warmup(arr);
        BlobBenchmark::benchmark(arr);
    }

    MPI_Finalize();
    return 0;
}
