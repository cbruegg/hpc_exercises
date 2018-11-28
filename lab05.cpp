#include <memory>
#include <vector>
#include <ostream>
#include <random>
#include <iostream>
#include <thread>
#include <fstream>
#include <chrono>
#include <pthread.h>

using namespace std;

class Main final {
public:
    static int main(const int argc, char *const argv[]) {
        if (argc < 2) {
            cerr << "Arguments: tosses threads" << endl;
            return EXIT_FAILURE;
        }

        const auto tosses = stou64(argv[1]);
        const auto totalRanks = stoi(argv[2]);
        const auto tossesPerRank = tosses / totalRanks;

        chrono::steady_clock::time_point threadLocalBegin = std::chrono::steady_clock::now();
        uint64_t globalTossesInCircleResultThreadLocal = threadLocalVersion(totalRanks, tossesPerRank);
        chrono::steady_clock::time_point threadLocalEnd = std::chrono::steady_clock::now();

        chrono::steady_clock::time_point sharedBegin = std::chrono::steady_clock::now();
        uint64_t globalTossesInCircleResultShared = sharedVersion(totalRanks, tossesPerRank);
        chrono::steady_clock::time_point sharedEnd = std::chrono::steady_clock::now();

        const auto piEstimateThreadLocal = (4.0 * globalTossesInCircleResultThreadLocal) / tosses;
        cout << globalTossesInCircleResultThreadLocal
             << " tosses landed inside the circle in thread-local version. π ≈ "
             << piEstimateThreadLocal
             << endl;

        const auto piEstimateShared = (4.0 * globalTossesInCircleResultShared) / tosses;
        cout << globalTossesInCircleResultShared
             << " tosses landed inside the circle in shared version. π ≈ "
             << piEstimateShared
             << endl;

        cout << "Times: threadLocal shared" << endl;
        cout << chrono::duration_cast<chrono::milliseconds>(threadLocalEnd - threadLocalBegin).count()
             << " "
             << chrono::duration_cast<chrono::milliseconds>(sharedEnd - sharedBegin).count()
             << endl;

        return 0;
    }

private:

    static pthread_mutex_t mutex;
    static uint64_t numberOfTossesSoFar;
    static uint64_t tossesWithinCircle;

    static uint64_t sharedVersion(const int totalRanks, const uint64_t tossesPerRank) {
        pthread_mutex_init(&mutex, nullptr);
        const uint64_t totalTosses = tossesPerRank * totalRanks;

        auto *const threads = (pthread_t *) malloc(totalRanks * sizeof(pthread_t));
        for (unsigned i = 0; i < totalRanks; i++) {
            pthread_create(&threads[i], nullptr, sharedThreadFuncVersion, (void *) &totalTosses);
        }
        for (unsigned i = 0; i < totalRanks; i++) {
            pthread_join(threads[i], nullptr);
        }
        free(threads);

        return tossesWithinCircle;
    }

    static void *sharedThreadFuncVersion(void *const totalTossesPtr) {
        const auto totalTosses = *static_cast<uint64_t *>(totalTossesPtr);
        uniform_real_distribution<double> uniform(-1, +1);
        default_random_engine engine;

        while (true) {
            const auto x = uniform(engine);
            const auto y = uniform(engine);
            const auto distanceSquared = x * x + y * y;

            pthread_mutex_lock(&mutex);
            auto numberOfTossesSoFar = Main::numberOfTossesSoFar;
            auto tossesWithinCircle = Main::tossesWithinCircle;

            if (numberOfTossesSoFar < totalTosses) {
                if (distanceSquared <= 1) {
                    tossesWithinCircle++;
                }
                numberOfTossesSoFar++;

                Main::numberOfTossesSoFar = numberOfTossesSoFar;
                Main::tossesWithinCircle = tossesWithinCircle;

                pthread_mutex_unlock(&mutex);
            } else {
                pthread_mutex_unlock(&mutex);
                break;
            }
        }

        return nullptr;
    }

    static uint64_t threadLocalVersion(const int totalRanks, const uint64_t tossesPerRank) {
        auto *const threads = (pthread_t *) malloc(totalRanks * sizeof(pthread_t));
        for (unsigned i = 0; i < totalRanks; i++) {
            pthread_create(&threads[i], nullptr, threadLocalThreadFuncVersion, (void *) &tossesPerRank);
        }

        uint64_t globalTossesInCircleResult = 0;
        for (unsigned i = 0; i < totalRanks; i++) {
            void *retVal;
            pthread_join(threads[i], &retVal);

            globalTossesInCircleResult += *static_cast<uint64_t *>(retVal);
            free(retVal);
        }

        free(threads);
        return globalTossesInCircleResult;
    }

    static void *threadLocalThreadFuncVersion(void *const totalTossesPtr) {
        uniform_real_distribution<double> uniform(-1, +1);
        default_random_engine engine;

        const auto totalTosses = *static_cast<uint64_t *>(totalTossesPtr);
        uint64_t result = 0;
        for (uint64_t toss = 0; toss < totalTosses; toss++) {
            const auto x = uniform(engine);
            const auto y = uniform(engine);
            const auto distanceSquared = x * x + y * y;
            if (distanceSquared <= 1) {
                result++;
            }
        }

        auto resultPtr = (uint64_t *) malloc(sizeof(uint64_t));
        *resultPtr = result;
        return resultPtr;
    }

    static uint64_t stou64(const string &s) {
        const unsigned long long result = stoull(s);
        if (result > numeric_limits<uint64_t>::max()) {
            throw invalid_argument("Number is too large!");
        }
        return result;
    }

    static uint64_t readTossesFromCinOrExit() {
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
};

pthread_mutex_t Main::mutex = {};
uint64_t Main::numberOfTossesSoFar = 0;
uint64_t Main::tossesWithinCircle = 0;

int main(int argc, char *argv[]) {
    return Main::main(argc, argv);
}