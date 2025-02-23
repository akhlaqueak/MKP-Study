all : 
	g++ -flto -Ofast -DNDEBUG -march=native -o UMKP Graph.cpp -w -fopenmp
clean:
	rm -rf UMKP
debug:
	# g++ -g -DNDEBUG -o UMKP Graph.cpp -fopenmp
	g++ -g  -o UMKP Graph.cpp -fopenmp -gdwarf-4 
prof:
	g++ -pg -DNDEBUG -march=native -o UMKP-pg Graph.cpp -w -fopenmp