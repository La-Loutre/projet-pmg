#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define _XOPEN_SOURCE 600

#include "display.h"

unsigned sand[DIM][DIM];

// vecteur de pixel renvoyé par compute
struct {
  float R, G, B;
} couleurs[DIM][DIM];

unsigned get (unsigned x, unsigned y)
{
  return sand[y][x];
}

#if CASE == 1
// on met du sable dans chaque case
static void sand_init ()
{
  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      sand[y][x] = MAX_HEIGHT + 1;
    }
}
#endif // CASE == 0

#if CASE == 2
// on construit un seul gros tas de sable
static void sand_init ()
{
  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      sand[y][x] = 0;
    }
  sand[DIM/2][DIM/2] = 100000;
}
#endif // CASE == 1

void afficher()
{
  // NOTE: we don't print the edges
  for(int i = 1; i < DIM-1; i++) {
    for(int j = 1; j < DIM-1; j++) {
      printf("%2d ", sand[i][j]);
    }
    printf("\n");
  }
}

float *compute_eucl (unsigned iterations)
{
  for (unsigned i = 0; i < iterations; i++)
    {
      for (int y = 1; y < DIM-1; y++)
	{
	  for (int x = 1; x < DIM-1; x++)
	    if(sand[y][x] >= 4) {
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
    }
  return DYNAMIC_COLORING;
}

float *compute_naive (unsigned iterations)
{
  for (unsigned i = 0; i < iterations; i++)
    {
      for (int y = 1; y < DIM-1; y++)
	{
	  for (int x = 1; x < DIM-1; x++)
	    if(sand[y][x] >= 4) {
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
    }
  return DYNAMIC_COLORING;
}

float *compute_omp (unsigned iterations)
{
  for (unsigned i = 0; i < iterations; i++)
    {
#pragma omp parallel for collapse(2) num_threads(4)
      for (int y = 1; y < DIM-1; y++)
	{
	  for (int x = 1; x < DIM-1; x++)
	    if(sand[y][x] >= 4) {
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
    }
  return DYNAMIC_COLORING;
}


int main (int argc, char **argv)
{
  printf("DIM %d CASE %d\n", DIM, CASE);

  sand_init ();

#if METHOD == 1
  display_init (argc, argv,
		DIM,              // dimension ( = x = y) du tas
		MAX_HEIGHT,       // hauteur maximale du tas
		get,              // callback func
		compute_naive);   // callback func
#endif // METHOD seq
#if METHOD == 2
  display_init (argc, argv,
		DIM,              // dimension ( = x = y) du tas
		MAX_HEIGHT,       // hauteur maximale du tas
		get,              // callback func
		compute_omp);     // callback func
#endif // METHOD openmp

  afficher();

  return 0;
}
