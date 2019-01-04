#pragma once

#include <stdio.h>
#include <stdlib.h>

class Playfield;

class Car {
public:
    Car(int x_, int y_, int dir_, int len_);

    bool intersects(const Car &) const;

    bool onField(int szx, int szy) const;

    bool checkWon(int cond, int szx, int szy) const;

    Car moved(bool forward) const;

    size_t hash() const {
        return size_t(x ^ (y
                << 1)); //Only those two things need be checked, the rest never changes (and is thus a bad criterion to hash by)

    };

    bool operator==(const Car &other) const {
        if (x == other.x &&
            y == other.y) { return true; } //Only those two things need be checked, the rest never changes!
        return false;
    };

    int getDir() const { return dir; };

private:
    int x;//Lower left point
    int y;
    int dir; // 0 hor, 1 vert
    int len;
};