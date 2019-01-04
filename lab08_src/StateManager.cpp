#include "StateManager.h"

#include <mutex>
#include <algorithm>

using namespace std;

void StateManager::enterSolution(const State &res) {
    if (res.solutionSize() < best_solution_size) {
        printf("Solution updated, size %d \n", res.solutionSize());
        best_solution = res;
        best_solution_size = res.solutionSize();
    }
}

bool StateManager::claim(const State &val) {


    /*
        If we find that some other task has already claimed the exact configuration and reached it in as many or less steps,
        we return false so it won't be worked on anymore.
    */
    auto val_it = state_set.find(val);
    if (val_it != state_set.cend() && val.solutionSize() >= val_it->solutionSize()) {
        return false; //Someone else is already doing that / has already done that
    }

    //Two val elements are treated as equal, even if the steps to reach them are different.
    //Erasing and re-inserting updates that value
    state_set.erase(val);
    state_set.insert(val);


    return true;
}


void StateManager::printBestSolution() {
    best_solution.printSolution();
}

