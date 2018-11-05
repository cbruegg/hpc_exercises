#include <mpi.h>

namespace hello {

    void sayHello() {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        std::cout << "Hello! I'm rank " << rank << " out of " << size << '\n';
    }

}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    hello::sayHello();
    MPI_Finalize();
    return 0;
}
