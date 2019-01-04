#include "Car.h"


#include "State.h"//for Playfield class


using namespace std;

Car::Car(int x_, int y_, int dir_, int len_)
        : x(x_), y(y_), dir(dir_), len(len_) {}

bool Car::intersects(const Car &other) const {
    for (int i = 0; i < len; i++) {
        int tx = x + i * (1 - dir);
        int ty = y + i * (dir);
        for (int oi = 0; oi < other.len; oi++) {
            int otx = other.x + oi * (1 - other.dir);
            int oty = other.y + oi * (other.dir);
            if (tx == otx && ty == oty) {
                //printf("Intersection at %d %d\n",tx,ty);
                return true;
            }
        }

    }
    return false;
}

//Returns a copy of this car that is identical except for having moved by 1 tile in the specified direction
Car Car::moved(bool forward) const {
    Car res(*this); //copy this car

    if (forward) {
        if (dir == 0) { ++res.x; } else { ++res.y; }
    } else {
        if (dir == 0) { --res.x; } else { --res.y; }
    }

    return res;
}


bool Car::onField(int szx, int szy) const {
    if (x < 0 || y < 0) { return false; }
    if (dir == 0 && x + len - 1 >= szx) { return false; }
    if (dir == 1 && y + len - 1 >= szy) { return false; }
    return true;
}


bool Car::checkWon(int cond, int szx, int szy) const {
    switch (cond) { //0 left, 1 top, 2 right 3 bot
        case 0:
            if (x == 0) return true;
            break;
        case 1:
            if (y + len == szy) return true;
            break;
        case 2:
            if (x + len == szx) return true;
            break;
        case 3:
            if (y == 0) return true;
            break;

    }
    return false;


}