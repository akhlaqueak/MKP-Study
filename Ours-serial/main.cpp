#include "Graph.h"
#include "Utility.h"
#include "Timer.h"

using namespace std;
#define LEN_LIMIT (1<<10)
char filename[LEN_LIMIT];

#define DENSE_SEARCH

int main(int argc, char *argv[]) {
	printf("\n-----------------------------------------------------------------------------------------\n");
	if (argc == 3) {
	#ifndef DENSE_SEARCH
		strncpy(filename, argv[1], LEN_LIMIT);
		int k = atoi(argv[2]);
		Graph *graph = new Graph(filename, k);
		graph->read();
		graph->search();
		graph->write();
		// delete graph; // there are some bugs in releasing memory
	#else
		strncpy(filename, argv[1], LEN_LIMIT);
		int k = atoi(argv[2]);
		Graph *graph = new Graph(filename, k);
		graph->read();
		graph->search();
		graph->search_dense();
		graph->write();
		// delete graph; // there are some bugs in releasing memory
	#endif
	}
	else printf("[usage]: exe file k [dense]\n");
	printf("-----------------------------------------------------------------------------------------\n\n");
}