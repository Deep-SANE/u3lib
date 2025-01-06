include config.mk

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
PRGS = $(basename $(OBJS))

#CXXFLAGS += -I../include

#LDFLAGS := -L../src -lSU3 $(LDFLAGS)

# debug messages:
#$(info $$SRCS is [$(SRCS)])
#$(info $$OBJS is [$(OBJS)])
#$(info $$PRGS is [$(PRGS)])
#$(info $$DEPS is [$(DEPS)])

.PHONY: all
all: $(PRGS)

$(PRGS): %: %.o
	$(CXX) -o $@ $< $(LDFLAGS)

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -rf $(OBJS) $(PRGS)
