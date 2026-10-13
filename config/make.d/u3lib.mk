# Include paths for u3lib's own headers and its vendored utils core.
# U3LIB_ROOT is set by config.mk (resolved from the tree layout).
CXXFLAGS += -I$(U3LIB_ROOT)/include
CXXFLAGS += -I$(U3LIB_ROOT)/include/internal
CXXFLAGS += -I$(U3LIB_ROOT)/utils
