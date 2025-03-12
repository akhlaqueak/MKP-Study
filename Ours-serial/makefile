all : 
	echo -e "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP main.cpp Graph.cpp -w
clean:
	rm -rf UMKP
debug:
	g++ -g -DNDEBUG -o UMKP main.cpp Graph.cpp
prof:
	g++ -pg -DNDEBUG -march=native -o UMKP-pg main.cpp Graph.cpp -w

versions:
	echo -e "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-Default main.cpp Graph.cpp -w

	echo -e "\n\
	#define _SECOND_ORDER_PRUNING_\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-ctcp main.cpp Graph.cpp -w

	echo -e "\n\
	#define SR_BRANCHING R_branching\n\
	#define PART_BRANCH (true)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-R_Br main.cpp Graph.cpp -w

	echo -e "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (true)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-S_Br main.cpp Graph.cpp -w

	echo -e "\n\
	#define SR_BRANCHING SR_branching\n\
	#define PART_BRANCH (true)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-SR_Br main.cpp Graph.cpp -w

	echo -e "\n\
	#define SR_BRANCHING SR_branching\n\
	// the SR_branching switch does not matter when PART_BRANCH is false\n\
	#define PART_BRANCH (false)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-Pivot_Br main.cpp Graph.cpp -w

	echo -e "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n\
	#define SR_BOUND" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-SR_Bound main.cpp Graph.cpp -w

	echo -e "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n\
	#define R_BOUND" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-R_Bound main.cpp Graph.cpp -w

	echo -e "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n\
	#define S_BOUND" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-S_Bound main.cpp Graph.cpp -w


	echo -e "\n\
	#define SR_BRANCHING S_branching\n\
	#define DIS_BR1\n\
	#define PART_BRANCH (K<=5)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-dis-BR1 main.cpp Graph.cpp -w


	echo -e "\n\
	#define SR_BRANCHING S_branching\n\
	#define DIS_BR2\n\
	#define PART_BRANCH (K<=5)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-dis-BR2 main.cpp Graph.cpp -w


	echo -e "\n\
	#define S2RULE\n\
	#define DIS_BR1\n\
	#define DIS_BR2\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-dis-BR12 main.cpp Graph.cpp -w


	echo -e "\n\
	#define S2RULE\n\
	#define DIS_UBR2\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-dis-UBR2 main.cpp Graph.cpp -w