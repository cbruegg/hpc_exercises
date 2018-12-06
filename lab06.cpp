#include<stdio.h>
#include<stdlib.h>
#include <string>
#include <iostream>
#include <memory>
#include <iomanip>

using namespace std;

const unsigned char LIVE = 1;
const unsigned char DEAD = 0;

void writePpm(const shared_ptr<unsigned char[]> gen, const string &filename, const size_t sideLength) {
    FILE *fp = fopen(filename.c_str(), "wb");
    fprintf(fp, "P6\n%zd %zd\n255\n", sideLength, sideLength);

    const unsigned char black[] = {0x00, 0x00, 0x00};
    const unsigned char white[] = {0xff, 0xff, 0xff};

    for (size_t i = 0; i < sideLength * sideLength; i++) {
        fwrite(gen[i] == 1 ? black : white, 1, 3, fp);
    }

    fclose(fp);
}


const shared_ptr<unsigned char[]> readPpm(const string &filename, size_t *const sideLength) {

    size_t h = 0;
    const auto fp = fopen(filename.c_str(), "rb");
    fscanf(fp, "P6\n%zu %zu\n255\n", sideLength, &h);

    if (*sideLength != h) exit(1);

    const shared_ptr<unsigned char[]> gen(new unsigned char[(*sideLength) * (*sideLength)]);

    for (auto i = 0ul; i < (*sideLength) * (*sideLength); i++) {
        unsigned char buf[3];
        fread(buf, 1, 3, fp);

        gen[i] = buf[0] == 0 && buf[1] == 0 && buf[2] == 0 ? LIVE : DEAD;
    }

    fclose(fp);

    return gen;
}

void update(const shared_ptr<unsigned char[]> data, size_t sideLength) {

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
    const auto prev = readPpm(fileName, &sideLength);
    for (auto i = 0; i < iterations; i++) {
        update(prev, sideLength);

        stringstream outFile;
        outFile << "output-" << setfill('0') << setw(2) << i << ".ppm";
        writePpm(prev, outFile.str(), sideLength);
    }
    return 0;
}