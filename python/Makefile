#BIN = ../bin

SRCS = $(wildcard *.cpp)
#PRGS = $(addprefix $(BIN)/,$(basename $(SRCS)))

.PHONY: all
all: $(SRCS)
	$(CXX) $(CXXFLAGS) -fPIC -shared -o libu3lib.so $(SRCS)

#.PHONY: clean
#clean:
#	rm -rf $(PRGS)
