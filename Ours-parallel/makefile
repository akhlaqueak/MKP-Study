UNAME_S := $(shell uname -s)

# Default flags
CXXFLAGS = -Ofast -DNDEBUG -march=native -w -std=c++17
LDFLAGS =

ifeq ($(UNAME_S),Darwin)
    # macOS (Apple Clang)
    CXX = clang++
    CXXFLAGS += -Xpreprocessor -fopenmp
    LDFLAGS += -lomp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib
else
    # Linux (GCC)
    CXX = g++
    CXXFLAGS += -fopenmp -flto
endif

all:
	$(CXX) $(CXXFLAGS) -o UMKP Graph.cpp $(LDFLAGS)

clean:
	rm -rf UMKP UMKP-pg

debug:
	$(CXX) -g -o UMKP Graph.cpp $(CXXFLAGS) $(LDFLAGS) -gdwarf-4

prof:
	$(CXX) -pg $(CXXFLAGS) -o UMKP-pg Graph.cpp $(LDFLAGS)