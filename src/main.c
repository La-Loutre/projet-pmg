#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#include "display.h"

#define _XOPEN_SOURCE 600

#define TEST 0
#define SEQEUCL 1
#define PARAOMP 2
#define FIVE_PILES 1
#define ONE_PILE 2

#define TIME_DIFF(t1, t2)						\
  ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

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
static void sand_init (sand_t sand)
{
  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      sand[y][x] = MAX_HEIGHT + 1;
    }
}
#endif

#if CASE == ONE_PILE
// on construit un seul gros tas de sable
static void sand_init (sand_t sand)
{
  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      sand[y][x] = 0;
    }
  sand[DIM/2][DIM/2] = 100000;
}
#endif

void afficher(sand_t sand)
{
  // NOTE: we don't print the edges
  for(int i = 1; i < DIM-1; i++) {
    for(int j = 1; j < DIM-1; j++) {
      printf("%2d ", sand[i][j]);
    }
    printf("\n");
  }
}

bool check(sand_t res, sand_t sand)
{
  for(int i = 0; i < DIM; i++) {
    for(int j = 0; j < DIM; j++) {
      if (res[i][j] != sand[i][j])
	return false;
    }
  }
  return true;
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


bool compute_eucl (sand_t sand)
{
  int changement = false;
  for (int y = 1; y < DIM-1; y++)
    {
      for (int x = 1; x < DIM-1; x++)
	if(sand[y][x] >= 4) {
	  changement = true;
	  int mod4 = sand[y][x] % 4;
	  int div4 = sand[y][x] / 4;
	  sand[y][x] = mod4;
	  if (y != 1)
	    sand[y-1][x] += div4;
	  if (y != DIM-2)
	    sand[y+1][x] += div4;
	  if (x != 1)
	    sand[y][x-1] += div4;
 	  if (x != DIM-2)
	    sand[y][x+1] += div4;
	}
    }
  //  return DYNAMIC_COLORING;
  return changement;
}

bool compute_naive (sand_t sand)
{
  bool changement = false;
  for (int y = 1; y < DIM-1; y++)
    {
      for (int x = 1; x < DIM-1; x++)
	if(sand[y][x] >= 4) {
	  changement = true;
	  sand[y][x] -= 4;
	  if (y != 1)
	    sand[y-1][x] += 1;
	  if (y != DIM-2)
	    sand[y+1][x] += 1;
	  if (x != 1)
	    sand[y][x-1] += 1;
	  if (x != DIM-2)
	    sand[y][x+1] += 1;
	}
    }
  return changement;
  //  return DYNAMIC_COLORING;
}

bool compute_omp (sand_t sand)
{
  bool changement = false;
#pragma omp parallel reduction(||:changement)
#pragma omp for collapse(2)
  for (int y = 1; y < DIM-1; y++)
    {
      for (int x = 1; x < DIM-1; x++)
	if(sand[y][x] >= 4) {
	  changement = true;
	  sand[y][x] -= 4;
	  if (y != 1)
	    sand[y-1][x] += 1;
	  if (y != DIM-2)
	    sand[y+1][x] += 1;
	  if (x != 1)
	    sand[y][x-1] += 1;
	  if (x != DIM-2)
	    sand[y][x+1] += 1;
	}
    }
  return changement;
  //  return DYNAMIC_COLORING;
}


int main (int argc, char **argv)
{
  omp_set_num_threads(4);
  printf("NTHREADS %d DIM %d CASE %d\n", omp_get_max_threads(), DIM, CASE);

  unsigned **sand = malloc(sizeof(unsigned*) * DIM);
  for(int i = 0; i < DIM; i++)
    sand[i] = malloc(DIM * sizeof(unsigned));
  unsigned **res = malloc(sizeof(unsigned*) * DIM);
  for(int i = 0; i < DIM; i++)
    res[i] = malloc(DIM * sizeof(unsigned));

  sand_init (sand);
  sand_init (res);

#if METHOD == SEQEUCL
  display_init (argc, argv,
		DIM,
		MAX_HEIGHT,
		get,
		iterate,
		compute_eucl,
		sand);
#endif // METHOD seq

#if METHOD == PARAOMP
  display_init (argc, argv,
		DIM,
		MAX_HEIGHT,
		get,
		iterate,
		compute_omp,
		sand);
#endif // METHOD openmp

#if METHOD == TEST
  struct timeval t1, t2;
  unsigned long seq_compute_time = 0;
  unsigned long compute_time = 0;

  gettimeofday (&t1, NULL);
  while(compute_naive(res));
  gettimeofday (&t2, NULL);
  seq_compute_time = TIME_DIFF(t1, t2);
  printf("SEQ %ld.%03ld ms\n", seq_compute_time/1000, seq_compute_time%1000);

  gettimeofday (&t1, NULL);
  while(compute_omp(sand));
  gettimeofday (&t2, NULL);
  printf("OMP %ld.%03ld ms\n", compute_time/1000, compute_time%1000);
  compute_time = 0;


  int err = 0;
  if (check(res, sand)) {
      fprintf(stderr,"OK OMP\n");
    }
  else {
      fprintf(stderr,"KO OMP\n");
      err = 1;
    }
  return err;
#endif
  return 0;
}
