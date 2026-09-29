# Build the u3lib Python extension (u3lib.so from pyu3lib.cpp).
#
# Dependency paths are NOT hardcoded here: they come from the qcx-config
# environment you source (qcx-config/env.d/<you>.env.sh), exactly like the C++
# build. Source that env, then from the u3lib root:  make python
# (or directly:  make -C python).
#
# The extension needs Boost/Eigen/wigxjpf + OpenMP, but NOT MPI/SLEPc/sns/libsam.
# The top-level Makefile exports the full app CXXFLAGS/LDFLAGS (which carry those),
# so we RESET them and re-include only the four fragments the extension needs.

CXXFLAGS := -std=c++20 -O2 -fPIC
LDFLAGS  :=

include $(QCX_CONFIG_DIR)/make.d/boost.mk      # -I$(BOOST_ROOT)/include
include $(QCX_CONFIG_DIR)/make.d/eigen.mk      # -I$(EIGEN_DIR)
include $(QCX_CONFIG_DIR)/make.d/wigxjpf.mk    # -DUSE_WIGXJPF, -I.../inc, -lwigxjpf
include $(QCX_CONFIG_DIR)/make.d/openmp.mk     # -DUSE_OPENMP, -fopenmp (libgomp)

# u3lib is header-only: its headers + the qcx-utils foundation.
CXXFLAGS += -I$(U3LIB_DIR)/include -I$(U3LIB_DIR)/include/internal
CXXFLAGS += -I$(QCX_UTILS_DIR)/include

# Build with plain g++ (the compiler mpicxx wraps): the extension does not use
# MPI, so we avoid linking libmpi into the .so. Override with `make CXX=...`.
ifdef OMPI_CXX
CXX := $(OMPI_CXX)
endif

# The interpreter you will `import u3lib` from (must have pybind11 installed).
# Override with `make PYTHON=/path/to/python`.
PYTHON ?= python3
PYBIND_INCLUDES := $(shell $(PYTHON) -m pybind11 --includes)

# Python extension link flags (macOS resolves Python symbols at load time).
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SOFLAGS := -shared -undefined dynamic_lookup
else
SOFLAGS := -shared
endif

# Header dependency tracking: g++ emits u3lib.d listing every header pyu3lib.cpp
# pulls in (u3lib + qcx-utils), so editing any of them triggers a rebuild.
DEPFLAGS := -MMD -MP -MF u3lib.d -MT u3lib.so

.PHONY: all clean
all: u3lib.so

u3lib.so: pyu3lib.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(PYBIND_INCLUDES) $(SOFLAGS) $< -o $@ $(LDFLAGS)

-include u3lib.d

clean:
	rm -f u3lib.so u3lib.d
