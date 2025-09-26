all : 	
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP Graph.cpp -w

CTCP : 	
	g++ -flto -Ofast -std=c++17 -DNDEBUG -march=native -o UMKP Graph.cpp -w -D_SECOND_ORDER_PRUNING_

clean:
	rm -rf UMKP
debug:
	g++ -std=c++17 -g -DNDEBUG -o UMKP main.cpp Graph.cpp
prof:
	g++ -std=c++17 -pg -DNDEBUG -march=native -o UMKP-pg main.cpp Graph.cpp -w
