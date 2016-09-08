CXX = g++ -std=c++11
SRC = lqr_solve.cpp

all: $(SRC)
	$(CXX) -o lqr_solve $^
