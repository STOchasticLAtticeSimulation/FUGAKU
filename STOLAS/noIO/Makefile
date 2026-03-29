MODEL = STOLAS

CXX := g++

CXXFLAGS := -std=c++26 -O3 -Xpreprocessor -fopenmp -march=native -ffast-math -flto -funroll-loops
CPPFLAGS := -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/boost/include

LDFLAGS := -L/opt/homebrew/opt/libomp/lib
LDLIBS := -lfftw3_threads -lfftw3 -lm -lomp

# LDFLAGS :=  -L/Users/theorygroup/.homebrew/opt/libomp/lib -L/Users/tomoakimurata/github_folder/usr/local/fftw/lib
# CPPFLAGS :=  -I/Users/theorygroup/.homebrew/opt/libomp/include -I/Users/tomoakimurata/github_folder/usr/local/boost/include -I/Users/tomoakimurata/github_folder/usr/local/fftw/include

all: $(MODEL)

$(MODEL): $(MODEL).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(MODEL).o: src/vec_op.hpp src/array_op.hpp STOLAS.hpp parameters.hpp model.hpp src/fft.hpp src/noise_bias.hpp src/zoom.hpp src/laplacian.hpp src/output.hpp


clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o
	$(RM) $(MODEL)
