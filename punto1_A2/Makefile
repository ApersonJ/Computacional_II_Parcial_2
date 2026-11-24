CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Iinclude

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:src/%.cpp=build/%.o)
HDR = $(wildcard include/*.h)

prog: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o prog

build/%.o: src/%.cpp $(HDR)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: prog
	./prog

clean:
	rm -rf build prog

.PHONY: run clean