#include<stdio.h>
#include<stdlib.h>
#include <string>
#include <iostream>
#include <memory>
#include <iomanip>
#include <vector>

using namespace std;

const unsigned char LIVE = 1;
const unsigned char DEAD = 0;

void writePpm(const vector<unsigned char> gen, const string &filename, const size_t sideLength) {
    FILE *fp = fopen(filename.c_str(), "wb");
    fprintf(fp, "P6\n%zd %zd\n255\n", sideLength, sideLength);

    const unsigned char black[] = {0x00, 0x00, 0x00};
    const unsigned char white[] = {0xff, 0xff, 0xff};

    for (size_t i = 0; i < sideLength * sideLength; i++) {
        fwrite(gen[i] == 1 ? black : white, 1, 3, fp);
    }

    fclose(fp);
}


const vector<unsigned char> readPpm(const string &filename, size_t *const sideLength) {

    size_t h = 0;
    const auto fp = fopen(filename.c_str(), "rb");
    fscanf(fp, "P6\n%zu %zu\n255\n", sideLength, &h);

    if (*sideLength != h) exit(1);

    vector<unsigned char> gen((*sideLength) * (*sideLength));

    for (auto i = 0ul; i < (*sideLength) * (*sideLength); i++) {
        unsigned char buf[3];
        fread(buf, 1, 3, fp);

        gen[i] = buf[0] == 0 && buf[1] == 0 && buf[2] == 0 ? LIVE : DEAD;
    }

    fclose(fp);

    return gen;
}

size_t countLiveNeighbors(const vector<unsigned char> data, size_t sideLength, size_t x, size_t y) {
    auto count = 0ul;
    for (auto neighborX = max(0ul, x - 1); neighborX <= min(sideLength - 1, x + 1); neighborX++) {
        for (auto neighborY = max(0ul, y - 1); neighborY <= min(sideLength - 1, y + 1); neighborY++) {
            if (neighborX == x && neighborY == y) continue;

            if (data[neighborY * sideLength + neighborX] == LIVE) {
                count++;
            }
        }
    }
    return count;
}

void update(const vector<unsigned char> &prev, vector<unsigned char> &next, size_t sideLength) {
    for (auto x = 0ul; x < sideLength; x++) {
        for (auto y = 0ul; y < sideLength; y++) {
            const auto liveNeighbors = countLiveNeighbors(prev, sideLength, x, y);

            const auto idx = y * sideLength + x;
            if (prev[idx] == LIVE) {
                if (liveNeighbors < 2 || liveNeighbors > 3) {
                    next[idx] = DEAD;
                }
            } else {
                if (liveNeighbors == 3) {
                    next[idx] = LIVE;
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "No initial population specified!" << endl;
        exit(1);
    } else if (argc < 3) {
        cerr << "No iteration count specified!" << endl;
        exit(1);
    }

    const string fileName = argv[1];
    const auto iterations = stoi(argv[2]);

    size_t sideLength;
    auto prev = readPpm(fileName, &sideLength);
    auto next = prev;

    for (auto i = 0; i < iterations; i++) {
        update(prev, next, sideLength);

        stringstream outFile;
        outFile << "output-" << setfill('0') << setw(2) << i << ".ppm";
        writePpm(next, outFile.str(), sideLength);

        prev = next;
    }
    return 0;
}