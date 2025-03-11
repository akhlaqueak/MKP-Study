all : 
	echo -e "\n\
	#define S2RULE\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5&&sparse)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP main.cpp Graph.cpp -w
clean:
	rm -rf UMKP
debug:
	g++ -g -DNDEBUG -o UMKP main.cpp Graph.cpp
prof:
	g++ -pg -DNDEBUG -march=native -o UMKP-pg main.cpp Graph.cpp -w

versions:
	echo -e "\n\
	#define S2RULE\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5&&sparse)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-Default main.cpp Graph.cpp -w

	echo -e "\n\
	#define S2RULE\n\
	#define SR_BRANCHING R_branching\n\
	#define PART_BRANCH (true)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-R_Br main.cpp Graph.cpp -w

	echo -e "\n\
	#define S2RULE\n\
	#define SR_BRANCHING SR_branching\n\
	#define PART_BRANCH (true)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-SR_Br main.cpp Graph.cpp -w

	echo -e "\n\
	#define S2RULE\n\
	#define SR_BRANCHING SR_branching\n\
	#define PART_BRANCH (false)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-Pivot_Br main.cpp Graph.cpp -w

	echo -e "\n\
	#define _SECOND_ORDER_PRUNING_\n\
	#define S2RULE\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5&&sparse)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-ctcp main.cpp Graph.cpp -w

	echo -e "\n\
	#define S2RULE\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5&&sparse)\n\
	#define SEESAW_BOUND" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-SR_Bound main.cpp Graph.cpp -w

	echo -e "\n\
	#define S2RULE\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5&&sparse)\n\
	#define COLOR_BOUND" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-R_Bound main.cpp Graph.cpp -w

	echo -e "\n\
	#define S2RULE\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5&&sparse)\n\
	#define PART_BOUND" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-S_Bound main.cpp Graph.cpp -w
