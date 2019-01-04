#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

#include <vector>

#include "Car.h"

class StateManager;


struct AMove {
    AMove(int a, int b) : what(a), how(b) {}

    int what;
    int how;
};


class State {
public:
    State()
            : cars(), steps(-1), movestack() {}

    State(std::vector <Car> cars_)
            : cars(cars_), steps(0), movestack() {};


    //Checks if this state is legal (no cars overlap, no car outside of the playing field)
    bool legal(const StateManager *pf) const;

    //Creates a followup state in which car nr. x is moved one tile forward/backward along its respective orientation
    State move_car(int car, bool forward) const;

    //How many cars are there?
    int carCount() const { return cars.size(); }

    //Checks if car 0 is touching the side of the playing field it is supposed to reach.
    bool won(const StateManager *pf) const;

    //Gives the number of steps it took to reach this configuration
    int solutionSize() const { return (int) movestack.size(); };


    //Checks if two States have the same confiuration of cars. The steps taken to reach them may be different.
    bool operator==(const State &other) const {
        for (unsigned int i = 0; i < cars.size(); ++i) {
            if (!(cars[i] == other.cars[i])) { return false; }
        }
        return true;
    };

    //Print the steps that brought us to this configuration
    void printSolution() const;

    //For std::unordered_map
    size_t hash() const;

private:
    std::vector <Car> cars;
    int steps;
    std::vector <AMove> movestack;
};


struct StateHash {
    std::size_t operator()(const State &k) const {
        return k.hash();
    }
};

struct StateEqual {
    bool operator()(const State &lhs, const State &rhs) const {
        return lhs == rhs;
    }
};



