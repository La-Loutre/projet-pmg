#ifndef _DISPLAY_H_IS_DEF_
#define _DISPLAY_H_IS_DEF_

#include <stdbool.h>

typedef unsigned **sand_t;

typedef unsigned (*get_func_t) (unsigned x, unsigned y,
				sand_t sand);

typedef bool (*compute_func_t) (sand_t sand);

typedef float * (*iterate_func_t) (compute_func_t compute_func,
				   unsigned iterations,
				   sand_t sand);

#define STATIC_COLORING    ((float *)NULL)
#define DYNAMIC_COLORING   ((float *)1)

void display_init (int argc, char **argv,
		   unsigned dim, unsigned max_height,
		   get_func_t get_func,
		   iterate_func_t iterate_func,
		   compute_func_t compute_func,
		   sand_t sand);

#endif
