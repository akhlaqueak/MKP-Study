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
	#define SR_BRANCHING S_br\n\
	#define PART_BRANCH (K<=5&&sparse)" > switches.h
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP-R_br main.cpp Graph.cpp -w


