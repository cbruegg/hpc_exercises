#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <memory>
#include <vector>

using namespace std;

/*
	A block struct contains one "task" to do.
	It represents a square area of the fractal to be rendered with the upper left point at x,y and lengths of size.
	it also contains information as to where to insert the calculated pixels as a starting shift (targetPos) and a size in pixels (targetSize)
*/
class Block {
public:
    const double x;//Upper left point
    const double y;//Upper left point
    const double size;
    //Output pixels
    const int targetPos;
    const int targetSize;//Assumes a square!

    Block(double x, double y, double size, int targetPos, int targetSize) : x(x), y(y), size(size),
                                                                            targetPos(targetPos),
                                                                            targetSize(targetSize) {}
};


int checkMandelbrot(double real, double imag, int cutoff) {
    /*
       Task 2
       ------
       Implement the code for determining whether the point c (in the complex plane, i.e. c = real + i*imag) is in the Mandelbrot
               set (https://en.wikipedia.org/wiki/Mandelbrot_set).

       The point (c = real + i*imag) is in the Mandelbrot set if the following sequence of complex numbers z_n is bounded
               (i.e. |z_n| <= 2  which means  sqrt(Real(z_n)^2 + Imag(z_n)^2) <=2 )

                z_0 = 0
                z_1 = z_0^2 + c
                z_2 = z_1^2 + c

                ...

        i.e. z_{n+1} = (z_n)^2 + c

               Consider that z_n are complex numbers!


       Perform the iteration up to as many times as the cutoff states, at which point you assume the point is in the set.
       If the absolute value (|z_n|) is ever greater than 2, the series will diverge for sure and you can conclude that the point is not in the set.

               This function returns the number of iterations n after that sqrt(Re(z_n)^2 + Im(z_n)^2) > 2.


       Notes:
           - You can do this test without calling the expensive square root function.

       Solution:
           - Looks like the set in Wikipedia.
   */


    if (real > imag) return cutoff; //<-- This needs to be modified.
    return 0;                      //<-- This needs to be modified.

}


void HandleBlock(int myRank, Block block, MPI_Win const &window, int totalSizeX, vector<int> localResults,
                 int maxNumberIterations) {
    for (auto v = 0; v < block.targetSize; ++v) {
        for (auto b = 0; b < block.targetSize; ++b) {

            auto result = checkMandelbrot(
                    block.x + block.size * b / block.targetSize, block.y + block.size * v / block.targetSize,
                    maxNumberIterations
            );

            localResults[b + v * block.targetSize] = result;
        }
    }

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, window);
    for (auto v = 0; v < block.targetSize; ++v) {
        MPI_Put(localResults.data() + v * totalSizeX,
                totalSizeX,
                MPI_INT,
                0,
                block.targetPos + v * totalSizeX,
                totalSizeX,
                MPI_INT,
                window);

        /*
           Task 1
           ------
           Transfer the pixels calculated above to the memory exposed by Process 0 using MPI_Put

           Notes:
               - Due to the different layouts, you need to do this line by line
               - Target location for each line is the base position(block.targetPos) + an offset(v*totalSizeX) that depends on the dimensions of the target matrix
               - Lock and unlock. This causes the write to happen before we overwrite our local results with the next block. (Passive Target Synchronization).
           Solution:
               - A triangular shape appears as the output
       */
    }
    MPI_Win_flush_local(myRank, window);
    MPI_Win_unlock(0, window);
}

class Main {
public:
    static const constexpr unsigned char white[3] = {255, 255, 255};
    static const constexpr unsigned char black[3] = {0, 0, 0};
};


int main(int argc, char *argv[]) {
    /*
        We render a square area with side lengths size (3rd command line parameter) with its upper left point given by the first two parameters.
        A fourth parameter sets the maximum number of iterations, to possibly bring out finer details.
    */
    auto outputSizePixels = 2000;
    auto posX = -1.80;
    auto posY = -1.0;
    auto size = 2.0;
    auto maxNumberIterations = 100;

    // TODO remove
//    maxNumberIterations = 1;

    //MPI seems to do fine with providing these to all participating processes
    if (argc >= 4) {
        posX = strtod(argv[1], nullptr);
        posY = strtod(argv[2], nullptr);
        size = strtod(argv[3], nullptr);
    } else if (argc >= 5) {
        maxNumberIterations = stoi(argv[4]);
    }


    MPI_Init(&argc, &argv);
    int myRank, totalRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    MPI_Win window;
    int *shared_data = nullptr;

    /*
            Task 3
            ------
            Create and distribute multiple Blocks in such a way, that more than one process can work at a time and the likelihood of load imbalances is low

            Notes:
                - Each of the n processes should work on n blocks, which ideally are not oriented as a row or column to reduce the likelihood of load imbalances

            Solution:
                - Should yield a speedup over the previous version, scaling with the number of processes used.
                - Try higher numbers of iterations.
    */


    if (myRank == 0) {
        //Allocate memory and make available to others

        MPI_Aint siz = outputSizePixels * outputSizePixels * sizeof(int);


        MPI_Alloc_mem(siz, MPI_INFO_NULL, &shared_data);
        MPI_Win_create(shared_data, siz, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
    } else {
        //All other ranks make no memory available for remote access
        MPI_Win_create(shared_data, 0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
    }

    auto blocksPerDirection = 1;
//    auto totalBlocks = blocksPerDirection * blocksPerDirection;
    auto blockPixelSize = outputSizePixels / blocksPerDirection;
    auto bufferSize = blockPixelSize * blockPixelSize;

    {
        //Allocate a block of memory to keep our local results until sending
        auto myResultVals = vector<int>(bufferSize);

        if (myRank == 1) {
            auto blockCoordX = 0; //In full blocks
            auto blockCoordY = 0;

            auto x = posX + size / blocksPerDirection * blockCoordX; //upper left corner of the block we want to manage
            auto y = posY + size / blocksPerDirection * blockCoordY;
            auto blockSize = size / blocksPerDirection;

            auto targetSize = blockPixelSize;  //in pixels
            //Start of the first row in our target matrix.
            //Second row will start with an offset of +outputSizePixels, etc.
            auto targetPos = targetSize * (blockCoordX + blockCoordY * outputSizePixels);

            Block block(x, y, blockSize, targetPos, targetSize);
            HandleBlock(myRank, block, window, outputSizePixels, myResultVals, maxNumberIterations);
        }
    }


    //Before outputting the result we wait for all the values
    MPI_Barrier(MPI_COMM_WORLD);

    if (myRank == 0) {
        char filename[100];
        sprintf(filename, "Mandelbrot_x%f y%f size %f.ppm", posX, posY, size);
        auto fp = fopen(filename, "wb"); /* b - binary mode */
        fprintf(fp, "P6\n%d %d\n255\n", outputSizePixels, outputSizePixels);

        for (auto i = 0; i < outputSizePixels; ++i) {
            for (auto b = 0; b < outputSizePixels; ++b) {
                auto res = shared_data[i * outputSizePixels + b];
                if (res > 0) {
                    auto greyscaleValue = static_cast<unsigned char>(250 - 200 * res / maxNumberIterations);
                    unsigned char greyscale[3] = {greyscaleValue, greyscaleValue, greyscaleValue};
                    fwrite(greyscale, 1, 3, fp);
                } else {
                    fwrite(Main::black, 1, 3, fp);
                }

            }
        }

        fclose(fp);
    }


    MPI_Win_free(&window);
    MPI_Finalize();

    return 0;
}