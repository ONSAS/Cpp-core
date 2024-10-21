
CXX = g++

CXXFLAGS = -Wall -g -O2 -larmadillo

EXE = timeStepIteration.lnx

main: main.o solution.o
	$(CXX) -o $(EXE) main.o solution.o $(CXXFLAGS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp solution.cpp

solution.o: solution.cpp
	$(CXX) $(CXXFLAGS) -c solution.cpp

clean:
	rm $(EXE) *.o

# option direct from console without make:
#  g++ timeStepIteration.cpp -o timeStepIteration.lnx -O2 -larmadillo

