all : 	
	printf "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K <= 5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP Graph.cpp -w

clean:
	rm -rf UMKP
debug:
	g++ -std=c++17 -g -DNDEBUG -o UMKP main.cpp Graph.cpp
prof:
	g++ -std=c++17 -pg -DNDEBUG -march=native -o UMKP-pg main.cpp Graph.cpp -w

versions:
	printf "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-Default main.cpp Graph.cpp -w

	printf "\n\
	#define _SECOND_ORDER_PRUNING_\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-ctcp main.cpp Graph.cpp -w

	printf "\n\
	#define SR_BRANCHING R_branching\n\
	#define PART_BRANCH (true)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-R_Br main.cpp Graph.cpp -w

	printf "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (true)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-S_Br main.cpp Graph.cpp -w

	printf "\n\
	#define SR_BRANCHING SR_branching\n\
	#define PART_BRANCH (true)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-SR_Br main.cpp Graph.cpp -w

	printf "\n\
	#define SR_BRANCHING SR_branching\n\
	// the SR_branching switch does not matter when PART_BRANCH is false\n\
	#define PART_BRANCH (false)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-Pivot_Br main.cpp Graph.cpp -w

	printf "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n\
	#define SR_BOUND\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-SR_Bound main.cpp Graph.cpp -w

	printf "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n\
	#define R_BOUND\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-R_Bound main.cpp Graph.cpp -w

	printf "\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n\
	#define S_BOUND\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-S_Bound main.cpp Graph.cpp -w


	printf "\n\
	#define SR_BRANCHING S_branching\n\
	#define DIS_BR1\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-dis-BR1 main.cpp Graph.cpp -w


	printf "\n\
	#define SR_BRANCHING S_branching\n\
	#define DIS_BR2\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-dis-BR2 main.cpp Graph.cpp -w


	printf "\n\
	#define DIS_BR1\n\
	#define DIS_BR2\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-dis-BR12 main.cpp Graph.cpp -w

	printf "\n\
	#define DIS_UBR2\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-dis-UBR2 main.cpp Graph.cpp -w

	printf "\n\
	#define DIS_RR1\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-dis-RR1 main.cpp Graph.cpp -w

	printf "\n\
	#define DIS_RR2\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-dis-RR2 main.cpp Graph.cpp -w

	printf "\n\
	#define DIS_RR3\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-dis-RR3 main.cpp Graph.cpp -w

	printf "\n\
	#define DIS_CTCP\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-dis-CTCP main.cpp Graph.cpp -w


	printf "\n\
	#define BINARY_BRANCHING\n\
	#define SR_BRANCHING S_branching\n\
	#define PART_BRANCH (K<=5)\n" > switches.h
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP-Binary_Br main.cpp Graph.cpp -w

