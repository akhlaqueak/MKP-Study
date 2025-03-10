all : 
	echo "
	#define S2RULE
	#define PART_BRANCH (K<=5&&sparse)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP main.cpp Graph.cpp -w
clean:
	rm -rf UMKP
debug:
	g++ -g -DNDEBUG -o UMKP main.cpp Graph.cpp
prof:
	g++ -pg -DNDEBUG -march=native -o UMKP-pg main.cpp Graph.cpp -w

versions:


// #define _SECOND_ORDER_PRUNING_

// pruning switches
#define S2RULE

// SR_BRANCHING can take values S_branching, R_branching, SR_branching
#define SR_BRANCHING S_branching
// if PART_BRANCH is false, then pivot branch gets executed... 
// #define PART_BRANCH (false)
#define PART_BRANCH (K<=5&&sparse)


// Upper bounding switches... 
// #define SEESAW
// #define COLORBOUND
// #define PART_BOUND