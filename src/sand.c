#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include <sys/time.h>
#include <tgmath.h>
#include <assert.h>

#include "sand.h"
#include "display.h"

#define TIME_DIFF(t1, t2)						\
  ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

#define MAX_HEIGHT 4

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define FIVE_PILES 1
#define ONE_PILE 2

// vecteur de pixel renvoyé par compute
struct {
  float R, G, B;
} couleurs[DIM][DIM];

unsigned get (unsigned x, unsigned y, sand_t sand)
{
  return sand[y][x];
}

#if CASE == FIVE_PILES
// on met du sable dans chaque case
void sand_init (sand_t sand)
{
  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      if (y == 0 || y == DIM-1 || x == 0 || x == DIM-1)
	sand[y][x] = 0;
      else
	sand[y][x] = MAX_HEIGHT + 1;
    }
}
#endif

#if CASE == ONE_PILE
// on construit un seul gros tas de sable
void sand_init (sand_t sand)
{
  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      sand[y][x] = 0;
    }
  sand[DIM/2][DIM/2] = 100000;
}
#endif

void print_matrix(sand_t sand, int size)
{
  // NOTE: we don't print the edges
  for(int i = 1; i < size-1; i++) {
    for(int j = 1; j < size-1; j++) {
      printf("%2d ", sand[i][j]);
    }
    printf("\n");
  }
}

bool check_matrix(sand_t ref, sand_t sand)
{
  // NOTE: we don't check the edges
  for(int i = 1; i < DIM-1; i++) {
    for(int j = 1; j < DIM-1; j++) {
      if (ref[i][j] != sand[i][j])
	return false;
    }
  }
  return true;
}

void timeandcheck(char *name,
		  unsigned long ref_time,
		  unsigned long compute_time,
		  sand_t ref,
		  sand_t sand)
{
  fprintf(stderr, "%s %ld.%03ld ms ",
	  name, compute_time/1000, compute_time%1000);

  double speedup;
  if (ref_time == 0)
    speedup = 1;
  else {
    speedup = ref_time;
    speedup /= compute_time;
  }
  fprintf(stderr, "%.1f× ", speedup);

  if (check_matrix(ref, sand))
    fprintf(stderr,"%sSUCCESS%s\n", ANSI_COLOR_GREEN, ANSI_COLOR_RESET);
  else
    fprintf(stderr,"%sFAILURE%s\n", ANSI_COLOR_RED, ANSI_COLOR_RESET);
}

unsigned long process(char *name,
		      sand_t ref,
		      sand_t sand,
		      compute_func_t compute,
		      unsigned long ref_time,
		      bool loop,
		      int repeat)
{
  struct timeval t1, t2;
  unsigned long compute_time = 0;
  for (int i = 0; i < repeat; i++) {
    sand_init (sand);
    if (loop) {
      gettimeofday (&t1, NULL);
      while(compute(sand));
      gettimeofday (&t2, NULL);
    }
    else {
      gettimeofday (&t1, NULL);
      compute(sand);
      gettimeofday (&t2, NULL);
    }
    compute_time += TIME_DIFF(t1, t2);
  }
  compute_time /= repeat;
  timeandcheck(name, ref_time, compute_time, ref, sand);
  return compute_time;
}

float *iterate(compute_func_t compute_func,
	       unsigned iterations,
	       sand_t sand)
{
  for (unsigned i = 0; i < iterations; i++) {
    if (compute_func(sand) == false)
      break;
  }
  return DYNAMIC_COLORING;
}

sand_t create_sand_array_naive(int size)
{
  unsigned **sand_array = malloc(sizeof(unsigned*) * size);
  for(int i = 0; i < size; i++)
    sand_array[i] = malloc(size * sizeof(unsigned));
  return sand_array;

}

sand_t create_sand_array(int size)
{
  unsigned *raw_sand_array = malloc(sizeof(unsigned*) * size * size);
  unsigned **two_dim_sand_array = malloc(sizeof(unsigned**) * size);
  for (int i = 0; i < size; ++i)
    two_dim_sand_array[i] = &raw_sand_array[i*size];
  return two_dim_sand_array;

}