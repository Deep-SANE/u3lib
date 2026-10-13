# u3lib — self-contained build. No environment variables required:
#   make config     resolve all dependencies, write config/config.generated.mk
#   make            build the C++ programs
#   make test       build and run the regression suite
#   make python     build the Python extension (make python PYTHON=/path/to/python)
#   make install    install headers to a prefix (PREFIX=..., else prompts)
#
# `make config` is interactive; for scripted/CI use pass flags through, e.g.
#   make config CONFIG_ARGS="--yes"

U3LIB_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
CONFIG_DIR := $(U3LIB_ROOT)/config
export U3LIB_ROOT CONFIG_DIR

# Targets that must work before (or without) a generated configuration.
NOCONFIG_GOALS := config reconfigure clean distclean help

ifeq ($(filter $(NOCONFIG_GOALS),$(MAKECMDGOALS)),)
include $(CONFIG_DIR)/config.mk
endif
export

.PHONY: all
all:
	$(MAKE) -C src

.PHONY: config reconfigure
config reconfigure:
	@bash $(CONFIG_DIR)/configure.sh $(CONFIG_ARGS)

.PHONY: test
test:
	$(MAKE) -C test

.PHONY: python
python:
	$(MAKE) -C python

.PHONY: install
install: all
	@prefix="$(PREFIX)"; \
	if [ -z "$$prefix" ]; then \
	  read -r -p "Install prefix [$$HOME/.local]: " prefix </dev/tty || prefix=""; \
	  prefix="$${prefix:-$$HOME/.local}"; \
	fi; \
	echo "Installing u3lib headers to $$prefix/include"; \
	mkdir -p "$$prefix/include"; \
	cp -R $(U3LIB_ROOT)/include/. "$$prefix/include/"; \
	cp -R $(U3LIB_ROOT)/utils/. "$$prefix/include/"; \
	echo "Done. Downstream builds: add -I$$prefix/include"

# Clean removes build artifacts directly, without parsing the sub-Makefiles
# (those include dependency fragments that require a loaded configuration).
.PHONY: clean
clean:
	rm -rf $(U3LIB_ROOT)/bin
	rm -f $(U3LIB_ROOT)/python/u3lib.so $(U3LIB_ROOT)/python/u3lib.d

# distclean resets to a pristine checkout: removes the generated config and the
# external-libs/ tree that `make config` downloads and builds (it is not part of
# the repository).
.PHONY: distclean
distclean: clean
	rm -f $(CONFIG_DIR)/config.generated.mk
	rm -rf $(U3LIB_ROOT)/external-libs

# Build the manual reference-generation tool without running it.
.PHONY: reference-tool
reference-tool:
	$(MAKE) -C test reference-tool
