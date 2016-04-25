MAKE=make -C src
BINDIR=bin
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

test:
	$(MAKE) clean
	$(MAKE) METHOD=3 CASE=1 DIM=128
	$(MAKE) METHOD=3 CASE=1 DIM=512
	$(MAKE) METHOD=3 CASE=2 DIM=128
	$(MAKE) METHOD=3 CASE=2 DIM=512
	./bin/sandpiles-m3-c1-h4-d128
	./bin/sandpiles-m3-c2-h4-d128
	./bin/sandpiles-m3-c1-h4-d512
	./bin/sandpiles-m3-c2-h4-d512

clean:
	$(MAKE) clean

clear:
	$(MAKE) clear
