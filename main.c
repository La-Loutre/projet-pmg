#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _XOPEN_SOURCE 600

#include "display.h"

#ifndef DIM
#define DIM 128
#endif // DIM

#ifndef MAX_HEIGHT
#define MAX_HEIGHT 4
#endif // MAX_HEIGHT

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
static void sand_init ()
{
  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      sand[y][x] = MAX_HEIGHT + 1;
    }
}
#endif // CASE == 0

#if CASE == 2
static void sand_init ()
{
  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      sand[y][x] = 0;
    }
  sand[DIM/2][DIM/2] = 10000;
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


int main (int argc, char **argv)
{
  sand_init ();

  display_init (argc, argv,
		DIM,              // dimension ( = x = y) du tas
		MAX_HEIGHT,       // hauteur maximale du tas
		get,              // callback func
		compute_eucl);    // callback func
  //  afficher();
  return 0;
}
