# set if needed:
#EIGEN_DIR=/usr/include/eigen3
#WIGXJPF_DIR=$(HOME)/qcx/software/wigxjpf-1.13

CXX = g++

INC += -I$(WIGXJPF_DIR)/inc
#INC += -I$(BOOST_ROOT)/include
INC += -I$(EIGEN_DIR)

CXXFLAGS += -std=c++17
CXXFLAGS += $(INC)

# optimization / debugging options:
CXXFLAGS += -O3 
CXXFLAGS += -DNDEBUG 
#CXXFLAGS += -O0 
#CXXFLAGS += -g 
#CXXFLAGS += -fsanitize=address

LDFLAGS += -L$(WIGXJPF_DIR)/lib -lwigxjpf

## OpenMP if applicable:
#CXXFLAGS += -fopenmp
#LDFLAGS += -fopenmp
