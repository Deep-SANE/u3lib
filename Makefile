.PHONY: all
all: 
	$(MAKE) -C src

.PHONY: test
test:
	$(MAKE) -C src/test

.PHONY: clean
clean:
	$(MAKE) -C src clean
	$(MAKE) -C src/test clean
	
