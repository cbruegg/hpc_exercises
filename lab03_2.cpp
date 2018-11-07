#include <vector>
#include <memory>
#include <mpi.h>
#include <utility>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    cout << "sdfdf";

    MPI_Finalize();
    return 0;
}
