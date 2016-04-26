#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <tgmath.h>

#include "display.h"

#define _XOPEN_SOURCE 600

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define TEST 0
#define SEQEUCL 1
#define PAROMP 2
#define FIVE_PILES 1
#define ONE_PILE 2

#define MAX_HEIGHT 4

#define TIME_DIFF(t1, t2) \
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
      if (y == 0 || y == DIM-1 || x == 0 || x == DIM-1)
	sand[y][x] = 0;
      else
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

bool check(sand_t ref, sand_t sand)
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

  if (check(ref, sand))
    fprintf(stderr,"%sSUCCESS%s\n", ANSI_COLOR_GREEN, ANSI_COLOR_RESET);
  else
    fprintf(stderr,"%sFAILURE%s\n", ANSI_COLOR_RED, ANSI_COLOR_RESET);
}

unsigned long process(char *name,
		      sand_t ref,
		      sand_t sand,
		      compute_func_t compute,
		      unsigned long ref_time,
		      bool loop)
{
  struct timeval t1, t2;
  unsigned long compute_time = 0;
  sand_init (sand);
  gettimeofday (&t1, NULL);
  if (loop)
    while(compute(sand));
  else
    compute(sand);
  gettimeofday (&t2, NULL);
  compute_time = TIME_DIFF(t1, t2);
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

static int compute_eucl_chunk(sand_t sand)
{
  // We do expect memory to be continuous
  // doesn't work with naive sand
  unsigned *sand_one_array = &sand[0][0];
  bool changement = true;
  int chunk_stagned;
  int mod4;
  int div4;
  int block_size = 3;
  int nb_chunk = DIM / block_size;
  if (DIM % block_size != 0)
    nb_chunk += 1;

  int chunk[nb_chunk];
  memset(chunk, 1, nb_chunk);

  int start_offset;
  while(changement == true){
    changement = false;
    start_offset = DIM+1;
    for(int chunk_iter=0; chunk_iter < nb_chunk; ++chunk_iter)
      {
	if (chunk[chunk_iter]){
	  chunk_stagned = 0;
	  for (int cursor = DIM*block_size*chunk_iter+start_offset;
	       cursor < DIM*DIM && cursor < DIM*block_size*(chunk_iter+1)-1;
	       ++cursor)
	    {
	      if(sand_one_array[cursor] >= 4) {
		changement = true;
		chunk_stagned = 1;
		mod4 = sand_one_array[cursor] % 4;
		div4 = sand_one_array[cursor] / 4;
		sand_one_array[cursor] = mod4;
		sand_one_array[cursor - DIM] += div4;
		chunk[(cursor - DIM )/(DIM*block_size)] = 1;
		sand_one_array[cursor + DIM] += div4;
		chunk[(cursor + DIM )/(DIM*block_size)] = 1;
		sand_one_array[cursor - 1] += div4;
		sand_one_array[cursor + 1] += div4;
	      }
	    }

	  chunk[chunk_iter] = chunk_stagned;
	}
	start_offset = 0;
      }
  }
  return changement;
}

static inline int compute_eucl (sand_t sand)
{
  int changement = 0;
  int mod4;
  int div4;
  int chunk_size = 2;
  int nb_chunk = DIM / chunk_size ;
  int chunk[nb_chunk];
  for (int i=0;i<nb_chunk;++i)
    chunk[i]=1;
  int start=1;
  for (int chunk_iter = 0; chunk_iter < nb_chunk; ++chunk_iter){
    if (!chunk[chunk_iter])
      continue;
    chunk[chunk_iter] = 0;
    for (int y = chunk_iter*chunk_size+start; y < DIM-1 && y < chunk_iter*chunk_size+start+chunk_size; ++y)
      {
	for (int x = 1; x < DIM-1; ++x){
#if MAX_HEIGHT != 4
	if(sand[y][x] >= MAX_HEIGHT) {
	  changement = 1;
	  mod4 = sand[y][x] % MAX_HEIGHT;
	  div4 = sand[y][x] / MAX_HEIGHT;
	  sand[y][x] = mod4;
	  sand[y-1][x] += div4;
	  sand[y+1][x] += div4;
	  sand[y][x-1] += div4;
	  sand[y][x+1] += div4;
	}
#else
	switch(!(div4=sand[y][x] >> 2)){
	case 0:
	  changement = 1;
	  chunk[y/chunk_size] = 1;
	  chunk[(y-1)/chunk_size] |= div4;
	  chunk[(y+1)/chunk_size] |= div4;
	  sand[y][x] &= 3;
	  sand[y-1][x] += div4;
	  sand[y+1][x] += div4;
	  sand[y][x-1] += div4;
	  sand[y][x+1] += div4;
	}
#endif
      }
    }
    start=0;
  }
  return changement;
}

static inline int compute_naive (sand_t sand)
{
  int changement = 0;
  for (int y = 1; y < DIM-1; y++)
    {
      for (int x = 1; x < DIM-1; x++)
	if(sand[y][x] >= 4) {
	  changement = 1;
	  sand[y][x] -= 4;
	  sand[y-1][x] += 1;
	  sand[y+1][x] += 1;
	  sand[y][x-1] += 1;
	  sand[y][x+1] += 1;
	}
    }
  return changement;
}

static inline int compute_omp (sand_t sand)
{
  int changement = 0;
#pragma omp parallel shared(changement)
  {
    int nthreads = omp_get_num_threads();
    int myid = omp_get_thread_num();
    int chunk = DIM/nthreads;
    unsigned mysand [DIM][DIM]; // NOTE: base pointer inc, should be fast
    //memset(&mysand, 0, DIM*DIM);

    do {
#pragma omp barrier
#pragma omp single // barrier
      changement = 0;
#pragma omp for schedule(static, chunk) reduction(|:changement)
      for (int y = 1; y < DIM-1; y++) {
	for (int x = 1; x < DIM-1; x++) {
	  int val = sand[y][x];
#if MAX_HEIGHT != 4
	  if (val >= MAX_HEIGHT)
	    changement = 1;
#else
	  // NOTE: works only if MAX_HEIGHT == 4
	  changement = changement | (val >> 2);
#endif
	  val %= MAX_HEIGHT;
	  val += sand[y-1][x] / 4
	    + sand[y+1][x] / 4
	    + sand[y][x-1] / 4
	    + sand[y][x+1] / 4;
	  mysand[y][x] = val;
	}
      } // END PARALLEL FOR
      if (changement) {
#pragma omp for schedule(static, chunk)
	//for (int y = myid*chunk+!myid; y < (myid+1)*chunk-!myid; y++) {
	for (int y = 1; y < DIM-1; y++) {
	  for (int x = 1; x < DIM-1; x++) {
	    sand[y][x] = mysand[y][x];
	  }
	} // END PARALLEL FOR
      }
    } while(changement);
  } // END PARALLEL
  return changement;
}

static unsigned **create_sand_array_naive(int size)
{
  unsigned **sand_array = malloc(sizeof(unsigned*) * size);
  for(int i = 0; i < size; i++)
    sand_array[i] = malloc(size * sizeof(unsigned));
  return sand_array;

}
static unsigned **create_sand_array(int size)
{
  unsigned *raw_sand_array = malloc(sizeof(unsigned*) * size * size);
  unsigned **two_dim_sand_array = malloc(sizeof(unsigned**) * size);
  for (int i = 0; i < size; ++i)
    two_dim_sand_array[i] = &raw_sand_array[i*size];
  return two_dim_sand_array;

}

int main (int argc, char **argv)
{
  omp_set_num_threads(8);
  printf("BINDING %d ", omp_get_proc_bind());
  printf("NTHREADS %d DIM %d CASE %d\n", omp_get_max_threads(), DIM, CASE);

  unsigned **sand = create_sand_array(DIM);
  sand_init (sand);

#if METHOD == SEQEUCL
  display_init (argc, argv,
		DIM,
		MAX_HEIGHT,
		get,
		iterate,
		compute_eucl,
		sand);
  return 0;
#endif // METHOD SEQ EUCL

#if METHOD == PAROMP
  display_init (argc, argv,
		DIM,
		MAX_HEIGHT,
		get,
		iterate,
		compute_omp,
		sand);
  return 0;
#endif // METHOD PAR OMP

#if METHOD == TEST
  unsigned **ref = create_sand_array(DIM);
  unsigned long ref_time = 0;

  // NOTE: We use naive compute time for reference
  ref_time = process("SEQ REF", ref, ref, compute_naive, ref_time, true);

  ref_time = fmin(ref_time,
		  process ("SEQ EUCL", ref, sand, compute_eucl, ref_time, true));

  //process ("SEQ EUCL CHUNK", ref, sand, compute_eucl_chunk, ref_time, false);

  // NOTE: We use best sequential time for reference
  process ("PAR OMP", ref, sand, compute_omp, ref_time, false);

  fprintf(stderr,"\n");
  return 0;
#endif // METHOD TEST
}
