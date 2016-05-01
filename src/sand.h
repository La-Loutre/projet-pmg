#ifndef SAND_H
#define SAND_H

#include <stdbool.h>

#include "type.h"
typedef int (*check_func_type) (void *,void *);
unsigned get (unsigned x, unsigned y, sand_t sand);
void sand_init (sand_t sand);
void print_matrix(sand_t sand, int size);
int check_matrix(sand_t ref, sand_t sand);
void timeandcheck(char *name, unsigned long ref_time,
		  unsigned long compute_time,
		  void* ref,
		  void* sand,
		  int (*check_func) (void*, void*));
unsigned long process(char *name,
			sand_t ref,
			sand_t sand,
			compute_func_t compute,
			unsigned long ref_time,
			bool loop,
			int repeat);
float *iterate(compute_func_t compute_func, unsigned iterations, sand_t sand);
sand_t create_sand_array_naive(int size);
sand_t create_sand_array(int size);

#endif
