
#define _XOPEN_SOURCE 600

#include "display.h"

#include <stdio.h>
#include <stdlib.h>

//////////////////////////////////////////////////////////////////////////
// Tas de sable "fake" (juste pour tester)

#define DIM 128
#define MAX_HEIGHT 4
#include <math.h>

unsigned ocean[DIM][DIM];

// vecteur de pixel renvoyé par compute
struct {
  float R, G, B;
} couleurs[DIM][DIM];

// callback
unsigned get (unsigned x, unsigned y)
{
  return ocean[y][x];
}

// Tas de sable initial
static void sable_init ()
{
  unsigned dmax2 = MAX_HEIGHT;

  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      ocean[y][x] = 0;
    }
  ocean[DIM/2][DIM/2] = 1000000;
}

void afficher()
{
  // NOTE: we don't print the edges
  for(int i = 1; i < DIM-1; i++) {
    for(int j = 1; j < DIM-1; j++) {
      printf("%2d ", ocean[i][j]);
    }
    printf("\n");
  }
}

/*
// callback
float *compute (unsigned iterations)
{
  static int step = 0;
  for (unsigned i = 0; i < iterations; i++)
    {
      step++;
      for (int y = 0; y < DIM; y++)
	{
	  int v =  MAX_HEIGHT * (1+sin( 4* (y+step) * 3.14/ DIM)) / 4;
	  for (int x = 0; x < DIM; x++)
	    ocean[y][x]  = v;
	}
    }
  return DYNAMIC_COLORING; // altitude-based coloring
  // return couleurs;
}
*/

float *compute (unsigned iterations)
{
  for (unsigned i = 0; i < iterations; i++)
    {
      for (int y = 1; y < DIM-1; y++)
	{
	  for (int x = 1; x < DIM-1; x++)
	    if(ocean[y][x] >= 4) {
	      int mod4 = ocean[y][x] % 4;
	      int div4 = ocean[y][x] / 4;
	      ocean[y][x] = mod4;
	      if (y != 1)
		ocean[y-1][x] += div4;
	      if (y != DIM-2)
	      	ocean[y+1][x] += div4;
	      if (x != 1)
	      	ocean[y][x-1] += div4;
	      if (x != DIM-2)
	      	ocean[y][x+1] += div4;
	    }
	}
    }
  return DYNAMIC_COLORING;
}


int main (int argc, char **argv)
{
  sable_init ();

  display_init (argc, argv,
		DIM,              // dimension ( = x = y) du tas
		MAX_HEIGHT,       // hauteur maximale du tas
		get,              // callback func
		compute);         // callback func
  //  afficher();
  return 0;
}
