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
	g++ -pg -DNDEBUG -march=native -o UMKP-pg main.cpp Graph.cpp -w