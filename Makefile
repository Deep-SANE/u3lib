include $(QCX_CONFIG_DIR)/config.mk

export

.PHONY: all
all: 
	$(MAKE) -C src

.PHONY: test
test:
	$(MAKE) -C test

.PHONY: python
python:
	$(MAKE) -C python

.PHONY: clean
clean:
	$(MAKE) -C src clean
	$(MAKE) -C src/test clean
	$(MAKE) -C test clean
	$(MAKE) -C python clean
	

# Build the manual reference-generation tool without running it.
.PHONY: reference-tool
reference-tool:
	$(MAKE) -C test reference-tool
