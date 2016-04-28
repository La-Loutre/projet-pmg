#ifndef SAND_H
#define SAND_H

#include <stdbool.h>

#include "type.h"

unsigned get (unsigned x, unsigned y, sand_t sand);
void sand_init (sand_t sand);
void print_matrix(sand_t sand, int size);
bool check(sand_t ref, sand_t sand);
void timeandcheck(char *name, unsigned long ref_time,
		  unsigned long compute_time,
		  sand_t ref,
		  sand_t sand);
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
