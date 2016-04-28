#ifndef TYPE_H
#define TYPE_H

typedef unsigned **sand_t;

typedef unsigned (*get_func_t) (unsigned x, unsigned y,
				sand_t sand);

typedef int (*compute_func_t) (sand_t sand);

typedef float * (*iterate_func_t) (compute_func_t compute_func,
				   unsigned iterations,
				   sand_t sand);

#endif
