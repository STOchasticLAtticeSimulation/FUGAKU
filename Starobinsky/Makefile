MODEL = Starobinsky
SOURCE = ../source

#CXX := g++-13
CXX := g++
#CXX := icpc

#CXXFLAGS := -std=c++17 -O2 -fopenmp
CXXFLAGS := -std=c++14 -O2 -Xpreprocessor -fopenmp -lomp -lfftw3_threads -lfftw3 -lm
#CXXFLAGS := -std=c++11 -O2 -qopenmp

LDFLAGS := -L/opt/homebrew/opt/libomp/lib
CPPFLAGS := -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/boost/include

$(MODEL): $(MODEL).o $(SOURCE)/STOLAS_odeint.o
	$(CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(SOURCE)/STOLAS_odeint.o

$(MODEL).o: $(SOURCE)/STOLAS.hpp $(SOURCE)/parameters.hpp
STOLAS_odeint.o: $(SOURCE)/STOLAS.hpp $(SOURCE)/parameters.hpp $(SOURCE)/vec_op.hpp

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o
	$(RM) $(MODEL)
