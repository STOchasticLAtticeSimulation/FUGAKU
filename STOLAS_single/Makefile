MODEL = STOLAS
NOISE = noisemap
BIAS = biasmap
FUNC = functions

CXX := g++

CXXFLAGS := -std=c++14 -O2 -Xpreprocessor -fopenmp
CPPFLAGS := -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/boost/include

LDFLAGS := -L/opt/homebrew/opt/libomp/lib
LDLIBS := -lfftw3_threads -lfftw3 -lm -lomp

all: $(NOISE) $(BIAS) $(FUNC) $(MODEL)
$(NOISE): $(NOISE).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(BIAS): $(BIAS).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(FUNC): $(FUNC).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(MODEL): $(MODEL).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)


$(NOISE).o: parameters.hpp src/bias_noise.hpp src/vec_op.hpp src/fft.hpp
$(BIAS).o: parameters.hpp src/bias_noise.hpp src/vec_op.hpp src/fft.hpp
$(FUNC).o: parameters.hpp src/bias_noise.hpp src/vec_op.hpp src/fft.hpp
$(MODEL).o: src/vec_op.hpp STOLAS.hpp parameters.hpp model.hpp


clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o
	$(RM) $(NOISE)
	$(RM) $(BIAS)
	$(RM) $(FUNC)
	$(RM) $(MODEL)
