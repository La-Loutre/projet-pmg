MAKE=make -C src

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

clean:
	$(MAKE) clean

clear:
	$(MAKE) clear
