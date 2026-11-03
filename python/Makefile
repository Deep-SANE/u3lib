# Build the u3lib Python extension (u3lib.so from pyu3lib.cpp).
#
# Dependency paths are NOT hardcoded here: they come from config/config.generated.mk
# (written by `make config`), exactly like the C++ build. Build from the u3lib
# root with:  make python   (or:  make python PYTHON=/path/to/python).
#
# The top-level Makefile exports the full CXXFLAGS/LDFLAGS; we RESET them and
# re-include only the fragments the extension needs: Boost, Eigen, OpenMP,
# libwignernj, and krans (static) for the null-space solve. krans's wrapper falls
# back to dense QR, giving identical results via the shared RREF basis fixing.

CXXFLAGS := -std=c++20 -O2 -fPIC
LDFLAGS  :=

include $(CONFIG_DIR)/make.d/boost.mk      # -I$(BOOST_ROOT)/include
include $(CONFIG_DIR)/make.d/eigen.mk      # -I$(EIGEN_DIR)
include $(CONFIG_DIR)/make.d/openmp.mk     # -DUSE_OPENMP, -fopenmp (libgomp)

# Wigner (SU(2)) backend: libwignernj (the sole backend).
include $(CONFIG_DIR)/make.d/wignernj.mk   # USE_WIGNERNJ -> -lwignernj (shared, @rpath)

# Null-space solver: krans when available (static), else dense QR.
ifeq ($(USE_KRANS),1)
NULLSPACE_SOLVER := krans
LDFLAGS += -L$(KRANS_DIR)/build                  # krans static archives (built in place)
else
NULLSPACE_SOLVER := qr
endif
include $(CONFIG_DIR)/make.d/krans.mk        # USE_KRANS -> -lkrans_eigen -lkrans_core (static)
include $(CONFIG_DIR)/make.d/nullspace.mk  # -DNULLSPACE_SOLVER_KRANS, or default QR
include $(CONFIG_DIR)/make.d/pretty.mk     # quiet one-line output (make V=1 for full)

# u3lib is header-only: its headers + the utils foundation.
CXXFLAGS += -I$(U3LIB_ROOT)/include -I$(U3LIB_ROOT)/include/internal
CXXFLAGS += -I$(U3LIB_ROOT)/utils

# Build with plain g++ (the compiler mpicxx wraps): the extension does not use
# MPI, so we avoid linking libmpi into the .so. Override with `make CXX=...`.
ifdef OMPI_CXX
CXX := $(OMPI_CXX)
endif

# The interpreter you will `import u3lib` from (must have pybind11 installed).
# Override with `make PYTHON=/path/to/python`.
PYTHON ?= $(shell \
  best=""; fb=""; \
  for p in "$$CONDA_PREFIX/bin/python3" python3 "$$HOME/miniforge3/bin/python3" \
           python3.13 python3.12 python3.11 python3.10; do \
    command -v "$$p" >/dev/null 2>&1 || continue; \
    "$$p" -c 'import pybind11' >/dev/null 2>&1 || continue; \
    [ -z "$$fb" ] && fb="$$p"; \
    if "$$p" -c 'import ipykernel' >/dev/null 2>&1; then best="$$p"; break; fi; \
  done; \
  echo "$${best:-$${fb:-python3}}")
$(info [python] building for PYTHON=$(PYTHON))
PYBIND_INCLUDES := $(shell $(PYTHON) -m pybind11 --includes)
ifeq ($(strip $(PYBIND_INCLUDES)),)
$(warning PYTHON=$(PYTHON) has no pybind11; pass 'make python PYTHON=/path/to/python' with an interpreter that has pybind11 (e.g. a conda env or venv).)
endif

# Python extension link flags (macOS resolves Python symbols at load time).
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SOFLAGS := -shared -undefined dynamic_lookup
else
SOFLAGS := -shared
endif

# Header dependency tracking: g++ emits u3lib.d listing every header pyu3lib.cpp
# pulls in (u3lib + utils), so editing any of them triggers a rebuild.
DEPFLAGS := -MMD -MP -MF u3lib.d -MT u3lib.so

.PHONY: all clean
all: u3lib.so

u3lib.so: pyu3lib.cpp
	$(call pretty_print,CXX,$@)
	$(Q)$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(PYBIND_INCLUDES) $(SOFLAGS) $< -o $@ $(LDFLAGS)

-include u3lib.d

clean:
	rm -f u3lib.so u3lib.d
