MAKE= make -C src
BINDIR= bin

all:
	$(MAKE) clean
	$(MAKE) METHOD=1 CASE=1 DIM=128
	$(MAKE) clear
	$(MAKE) METHOD=1 CASE=1 DIM=512
	$(MAKE) clear
	$(MAKE) METHOD=1 CASE=2 DIM=128
	$(MAKE) clear
	$(MAKE) METHOD=1 CASE=2 DIM=512
	$(MAKE) clear
	$(MAKE) METHOD=2 CASE=1 DIM=128
	$(MAKE) clear
	$(MAKE) METHOD=2 CASE=1 DIM=512
	$(MAKE) clear
	$(MAKE) METHOD=2 CASE=2 DIM=128
	$(MAKE) clear
	$(MAKE) METHOD=2 CASE=2 DIM=512
	$(MAKE) clear
	$(MAKE) METHOD=3 CASE=1 DIM=128
	$(MAKE) clear
	$(MAKE) METHOD=3 CASE=1 DIM=512
	$(MAKE) clear
	$(MAKE) METHOD=3 CASE=2 DIM=128
	$(MAKE) clear
	$(MAKE) METHOD=3 CASE=2 DIM=512

test:
	$(MAKE) clean
	$(MAKE) METHOD=0 CASE=1 DIM=128
	$(MAKE) clear
	$(MAKE) METHOD=0 CASE=1 DIM=256
	$(MAKE) clear
	$(MAKE) METHOD=0 CASE=1 DIM=512
	$(MAKE) clear
	$(MAKE) METHOD=0 CASE=2 DIM=128
	$(MAKE) clear
	$(MAKE) METHOD=0 CASE=2 DIM=256
	$(MAKE) clear
	$(MAKE) METHOD=0 CASE=1 DIM=6
	$(MAKE) clear
	$(MAKE) METHOD=0 CASE=2 DIM=6
	./bin/sandpiles-m0-c1-d6
	./bin/sandpiles-m0-c2-d6
	./bin/sandpiles-m0-c1-d128
	./bin/sandpiles-m0-c2-d128
	./bin/sandpiles-m0-c1-d256
	./bin/sandpiles-m0-c2-d256
	./bin/sandpiles-m0-c1-d512
	./bin/sandpiles-m0-c2-d512


clean:
	$(MAKE) clean
