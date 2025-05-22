#BIN = ../bin

SRCS = $(wildcard *.cpp)
#PRGS = $(addprefix $(BIN)/,$(basename $(SRCS)))

CXXFLAGS += -shared -fPIC
CXXFLAGS += $(shell python -m pybind11 --includes)

.PHONY: all
all: $(SRCS)
	$(CXX) $(CXXFLAGS) -o u3lib.so $(SRCS) $(LDFLAGS)

#.PHONY: clean
#clean:
#	rm -rf $(PRGS)
