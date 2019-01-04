#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unordered_set>

#include <vector>


#include "StateManager.h"
#include "Car.h"

using namespace std;


/*
	state is the arrangement we want to check,
	manager is a pointer to the StateManager that coordinates the global state by keeping track of results and already checked states.
*/
void Check(State state, StateManager *manager) {
/*
	Task 1
	------
	Implement the body of this function. It should perform the following operations:
		-Immediately return if the State to be analysed took more steps to reach than our current best solution.
		-Try to claim its state in the manager. If the state is already claimed by another task, immediately return.
		-Check if the state is a winning state, if so enter its solution into the manager and return.
		-Iterate over all the cars (state.carCount()), for each car create the two followup states
		 created by moving the respective car forward or backward
		 (state.move_car(...) returns such a followup state from a given car number and direction)
		-Check whether the followup states created are legal states. If so recursively call Check(...) on them.
*/

    return;
}


int main(int argc, char *argv[]) {
/*
	Task 0
	------
	Draw with pen and paper the initial configuration built here.
	Work out a shortest solution to compare your algorithm against.
*/

    /*state_manager holds global information about the state of the solver */
    Playfield pf(8, 8, Playfield::GoalType::Right); //sizeX(8),sizeY(8),goal(2)
    StateManager *state_manager = new StateManager(pf);

    /*
        State holds information about one configuration we analyze. The initial state will be initialized with the following cars:
        Car(x,y,orientation,length): x,y coordinates of the lower left tile of the car,
                                     orientation 0: horizontally aligned, drives left/right
                                     orientation 1: vertically aligned, drives down/up
                                     length: how long is this car?
        The first car is the one that needs to reach the side of the playing field stated above, which is only possible if it has the proper orientation!
        The origin x=0, y=0 is the position in the lower left corner of the playing field.
    */
    vector <Car> cars_ = {Car(0, 4, 0, 2), Car(2, 4, 1, 3), Car(3, 2, 1, 3), Car(0, 2, 0, 2), Car(2, 1, 0, 2)};

/*
	Task 2
	------
	-Add OpenMP functionality to the solver by enclosing the main function's Check() call in a
	 #parallel environment and a taskgroup. Spawn the Check(...) function as a task.
	 Remember that the first Check(...) call only needs to be called by one thread
	-In Check() spawn more tasks for the recursion in a reasonable way.
	 Use the "default(none)" tag and explicitly state what a task can access and how it is accessed
*/

    Check(State(cars_), state_manager);

    state_manager->printBestSolution();

    delete state_manager;
    return 0;

}
