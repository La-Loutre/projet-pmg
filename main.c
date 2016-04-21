
#define _XOPEN_SOURCE 600

#include "display.h"

#include <stdio.h>
#include <stdlib.h>

//////////////////////////////////////////////////////////////////////////
// Tas de sable "fake" (juste pour tester)

#define DIM 128
#define MAX_HEIGHT  128
#include <math.h>

unsigned sable[DIM][DIM];

// vecteur de pixel renvoyé par compute
struct {
  float R, G, B;
} couleurs[DIM][DIM];

// callback
unsigned get (unsigned x, unsigned y)
{
  return sable[y][x];
}

// Tas de sable initial
static void sable_init ()
{
  unsigned dmax2 = MAX_HEIGHT;

  for (int y = 0; y < DIM; y++)
    for (int x = 0; x < DIM; x++) {
      sable[y][x] = MAX_HEIGHT / 4;
    }
}

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
	    sable[y][x]  = v;
	}
    }
  return DYNAMIC_COLORING; // altitude-based coloring
  // return couleurs;
}


int main (int argc, char **argv)
{
  sable_init ();

  display_init (argc, argv,
		DIM,              // dimension ( = x = y) du tas
		MAX_HEIGHT,       // hauteur maximale du tas
		get,              // callback func
		compute);         // callback func

  return 0;
}
