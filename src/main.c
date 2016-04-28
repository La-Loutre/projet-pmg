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
#include <emmintrin.h>
#include <immintrin.h>
#include "display.h"

#define _XOPEN_SOURCE 600

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define TEST 0
#define SEQEUCL 1
#define PAROMP 2
#define PAROMPSEM 3

#define FIVE_PILES 1
#define ONE_PILE 2

#define MAX_HEIGHT 4

#define TIME_DIFF(t1, t2)						\
  ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))


static sand_t create_sand_array(int size)
{
  unsigned *raw_sand_array = malloc(sizeof(unsigned*) * size * size);
  unsigned **two_dim_sand_array = malloc(sizeof(unsigned**) * size);
  for (int i = 0; i < size; ++i)
    two_dim_sand_array[i] = &raw_sand_array[i*size];
  return two_dim_sand_array;

}

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
		      bool loop,
		      int repeat)
{
  struct timeval t1, t2;
  unsigned long compute_time = 0;
  for (int i = 0; i < repeat; i++) {
    sand_init (sand);
    gettimeofday (&t1, NULL);
    if (loop)
      while(compute(sand));
    else
      compute(sand);
    gettimeofday (&t2, NULL);
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

static inline int compute_eucl_chunk (sand_t sand)
{
  int change = 0;
  int mod4;
  int div4;
  int chunk_size = 2;
  int nb_chunk = DIM / chunk_size ;
  int chunk[nb_chunk];
  for (int i=0;i<nb_chunk;++i)
    chunk[i]=1;
  int start=1;
  for (int chunk_iter = 0; chunk_iter < nb_chunk; ++chunk_iter) {
    if (!chunk[chunk_iter])
      continue;
    chunk[chunk_iter] = 0;
    for (int y = chunk_iter*chunk_size+start; y < DIM-1 && y < chunk_iter*chunk_size+start+chunk_size; ++y)
      {
	for (int x = 1; x < DIM-1; ++x){
#if MAX_HEIGHT != 4
	  if(sand[y][x] >= MAX_HEIGHT) {
	    change = 1;
	    mod4 = sand[y][x] % MAX_HEIGHT;
	    div4 = sand[y][x] / MAX_HEIGHT;
	    sand[y][x] = mod4;
	    sand[y-1][x] += div4;
	    sand[y+1][x] += div4;
	    sand[y][x-1] += div4;
	    sand[y][x+1] += div4;
	  }
#else
	  switch((div4=sand[y][x] >> 2)){
	  case 0:
	    continue;
	  default:
	    change = 1;
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
  return change;
}

static inline int compute_eucl (sand_t sand)
{

  int change = 0;
  int mod4;
  int div4;
  __m128i div4_simd,raw_value_simd,mod_simd,previousline_simd,
    nextline_simd,new_value_simd,simd_div,simd_div_mask,
    simd_mod_mask;
  int value_mod = 3;
  int value_div = 0xC0000000; //32 bits value
  int value_div_mask[4] ={value_div, 
			  value_div,
			  value_div,
			  value_div};
  int value_mod_mask[4] = {value_mod,
			   value_mod,
			   value_mod,
			   value_mod};
  simd_div_mask = _mm_loadu_si128(&value_div_mask[0]);
  simd_mod_mask = _mm_loadu_si128(&value_mod_mask[0]);
  
  for (int y = 1; y < DIM-1; ++y) {
    for (int x = 1; x < DIM-1; x+=4) {
#if MAX_HEIGHT == 4
      if(sand[y][x] >= MAX_HEIGHT) {
	change = 1;
	mod4 = sand[y][x] % MAX_HEIGHT;
	div4 = sand[y][x] / MAX_HEIGHT;
	sand[y][x] = mod4;
	sand[y-1][x] += div4;
	sand[y+1][x] += div4;
	sand[y][x-1] += div4;
	sand[y][x+1] += div4;
      }
#else

      for (int i=0;i<4;++i){
	div4=sand[y][x+i] >> 2;
	change |=div4;
      }
      
      //change |= div4; 
      /* SSE part */
      /*
       * SSE doesn't provide integer 
       * by integer division.
       * We have to trick
       * div4=sand[y][x] >> 2; for 4 elements
       */
      //Copy 4 element into 128b vector (integer)
      raw_value_simd = _mm_loadu_si128(&sand[y][x]);
      //Shift right 2 bits. 
      div4_simd = _mm_srli_si128(raw_value_simd,2);
      // We have to delete overflow by
      // applying a mask (xor 11000...11000...11000...11000)
      div4_simd = _mm_xor_si128(div4_simd,simd_div_mask);
      
      /* 
       * Now we have to do mod
       * We just apply a 3 mask already prepared
       * and its our new value
       * sand[y][x] &= 3; for 4 elements
       */
      new_value_simd = _mm_and_si128(raw_value_simd,
				     simd_mod_mask);

      /*
       * Upper and lower lines are easy to 
       * do because we never share y,x couple 
       * with our new_value_simd
       * sand[y-1][x] += div4;
       * sand[y+1][x] += div4;
       */
      previousline_simd = _mm_loadu_si128(&sand[y-1][x]);
      nextline_simd = _mm_loadu_si128(&sand[y+1][x]);
      previousline_simd = _mm_add_epi64(previousline_simd,
					div4_simd);
      nextline_simd = _mm_add_epi64(nextline_simd,
				    div4_simd);

      /*
       * For element on same line its different:
       * [... , x , x1  , x2   , x3  , ...]
       * [x-1 , x , x+1 , x+2  , x+3 , x+4]
       * [x   , x1, x+x2, x1+x3, x2  , x3 ]
       *
       * So we must store back data from vector 
       * register to memory
       */
      //Store back previous line
      _mm_storeu_si128(&sand[y-1][x],previousline_simd);
      //Store back next line
      _mm_storeu_si128(&sand[y+1][x],previousline_simd);
      //store back current line (4 elements)
      _mm_storeu_si128(&sand[y][x],new_value_simd);

      int x1_tmp,x2_tmp,x3_tmp;

      //First and last element can be done
      //directly (x-1 -> x+4)
      sand[y][x-1] += sand[y][x];
      sand[y][x+4] += sand[y][x+3];

      //For others line we must
      // do a swap before
      x1_tmp = sand[y][x+1];
      
      sand[y][x+1] += sand[y][x] + sand[y][x+2];
      sand[y][x] += x1_tmp;

      x2_tmp = sand[y][x+2];
      sand[y][x+2] += x1_tmp+sand[y][x+3];

      sand[y][x+3] += x2_tmp;

      //div4=sand[y][x] >> 2; 
      //change |= div4; 
      /* sand[y][x] &= 3; */
      /* sand[y-1][x] += div4; */
      /* sand[y+1][x] += div4; */
      /* sand[y][x-1] += div4; */
      /* sand[y][x+1] += div4; */
    
#endif
    }
  }
  return change;
}

static inline int compute_naive (sand_t sand)
{
  int change = 0;
  for (int y = 1; y < DIM-1; y++)
    {
      for (int x = 1; x < DIM-1; x++)
	if(sand[y][x] >= 4) {
	  change = 1;
	  sand[y][x] -= 4;
	  sand[y-1][x] += 1;
	  sand[y+1][x] += 1;
	  sand[y][x-1] += 1;
	  sand[y][x+1] += 1;
	}
    }
  return change;
}

static inline int compute_omp (sand_t sand)
{
  int change = 0;
#pragma omp parallel shared(change)
  {
    int nthreads = omp_get_num_threads();
    int myid = omp_get_thread_num();
    int chunk = (DIM-2)/nthreads;
    unsigned mysand [DIM][DIM]; // NOTE: base pointer offset, should be fast

    do {
#pragma omp barrier
#pragma omp single // barrier
      change = 0;

#pragma omp for schedule(static, chunk) reduction(|:change)
      for (int y = 1; y < DIM-1; y++) {
	for (int x = 1; x < DIM-1; x++) {
	  int val = sand[y][x];
#if MAX_HEIGHT != 4
	  if (val >= MAX_HEIGHT)
	    change = 1;
#else
	  // NOTE: works only if MAX_HEIGHT == 4
	  change = change | (val >> 2);
#endif
	  val %= MAX_HEIGHT;
	  val += sand[y-1][x] / 4
	    + sand[y+1][x] / 4
	    + sand[y][x-1] / 4
	    + sand[y][x+1] / 4;
	  mysand[y][x] = val;
	}
      } // END PARALLEL FOR

#pragma omp barrier
      // SYNCHRONISATION
      if (change) {
#pragma omp for schedule(static, chunk)
	for (int y = 1; y < DIM-1; y++) {
	  for (int x = 1; x < DIM-1; x++) {
	    sand[y][x] = mysand[y][x];
	  }
	} // END PARALLEL FOR
      }
    } while(change);
  } // END PARALLEL
  return change;
}

static inline int compute_omp_tile (sand_t sand)
{
  int change = 0;
#pragma omp parallel shared(change)
  {
    int nthreads = omp_get_num_threads();
    int myid = omp_get_thread_num();
    int chunk = (DIM-2)/nthreads/2;
    unsigned mysand [DIM][DIM]; // NOTE: base pointer offset, should be fast

    do {
#pragma omp barrier
#pragma omp single // barrier
      change = 0;

#pragma omp for schedule(static, chunk) reduction(|:change) collapse(2)
      for (int y = 1; y < DIM-1; y++) {
	for (int x = 1; x < DIM-1; x++) {
	  int val = sand[y][x];
#if MAX_HEIGHT != 4
	  if (val >= MAX_HEIGHT)
	    change = 1;
#else
	  // NOTE: works only if MAX_HEIGHT == 4
	  change = change | (val >> 2);
#endif
	  val %= MAX_HEIGHT;
	  val += sand[y-1][x] / 4
	    + sand[y+1][x] / 4
	    + sand[y][x-1] / 4
	    + sand[y][x+1] / 4;
	  mysand[y][x] = val;
	}
      } // END PARALLEL FOR

#pragma omp barrier
      // SYNCHRONISATION
      if (change) {
#pragma omp for schedule(static, chunk) collapse(2)
	for (int y = 1; y < DIM-1; y++) {
	  for (int x = 1; x < DIM-1; x++) {
	    sand[y][x] = mysand[y][x];
	  }
	} // END PARALLEL FOR
      }
    } while(change);
  } // END PARALLEL
  return change;
}

static inline int compute_omp_sem (sand_t sand)
{
  int niterations;
  int change = 0;
  int nthreads = omp_get_max_threads();
  sand_t aux = create_sand_array(DIM);

  // We will read the edges, so they should be set to 0
  memset(*aux, 0, DIM*DIM*sizeof(unsigned));

  sand_t swap[2] = {sand, aux};

  sem_t *locks = malloc(sizeof(sem_t)*(nthreads-1));
  for (int i = 0; i < nthreads-1; i++)
    assert(sem_init(&locks[i], 0, 0) ==0);

#pragma omp parallel shared(change)
  {
    sand_t read_from, write_to;
    read_from = swap[0];
    write_to = swap[1];

    int myid = omp_get_thread_num();
    int chunk = (DIM-2) / nthreads;
    int read = 0;
    int write = 1;

    do {
#pragma omp barrier
#pragma omp single // barrier
      change = 0;

#pragma omp for schedule(static, chunk) reduction(|:change)
      for (int y = 1; y < DIM-1; y++) {
	int chunk_number = (y-1) / chunk;
	// if nthreads is not a multiple of DIM
	// NOTE: one incrorrect branch prediction at maximum
	if (chunk_number >= nthreads)
	  chunk_number = nthreads -1;
	int first = chunk_number * chunk + 1;
	int last;
	if (chunk_number == nthreads-1)
	  last = DIM-2;
	else
	  last = first + chunk-1;

	// WAIT
	// NOTE: two incorrect branch predictions at maximum
	if (y == last && last != DIM-2) {
	  assert(sem_wait(&locks[chunk_number]) == 0);
	}
	for (int x = 1; x < DIM-1; x++) {
	  int val = read_from[y][x];

	  // UPDATE
	  // NOTE: works only if MAX_HEIGHT == 4
	  change = change | (val >> 2);
	  val &= 3 ;
	  val += read_from[y-1][x] / 4
	    + read_from[y+1][x] / 4
	    + read_from[y][x-1] / 4
	    + read_from[y][x+1] / 4;

	  write_to[y][x] = val;
	}
	// POST
	// NOTE: two incorrect branch predictions at maximum
	if (y == first && first != 1) {
	  assert(nthreads-1 > chunk_number-1);
	  assert(chunk_number >= 0);
	  assert (sem_post(&locks[chunk_number-1]) == 0);
	}

      } // END PARALLEL FOR
      read = 1 - read;
      write = 1 - write;
      read_from = swap[read];
      write_to = swap[write];

    } while(change);
 } // END PARALLEL
  free(*aux);
  free(aux);
  free(locks);
  return change;
}


 static sand_t create_sand_array_naive(int size)
 {
   unsigned **sand_array = malloc(sizeof(unsigned*) * size);
   for(int i = 0; i < size; i++)
     sand_array[i] = malloc(size * sizeof(unsigned));
   return sand_array;

 }


 int main (int argc, char **argv)
 {
   omp_set_num_threads(8);
   omp_set_nested(1);
   printf("BINDING %d ", omp_get_proc_bind());
   printf("NTHREADS %d DIM %d CASE %d\n", omp_get_max_threads(), DIM, CASE);

   sand_t sand = create_sand_array(DIM);
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

#if METHOD == PAROMPSEM
   display_init (argc, argv,
		 DIM,
		 MAX_HEIGHT,
		 get,
		 iterate,
		 compute_omp_sem,
		 sand);
   return 0;
#endif // METHOD PAR OMP SEM

#if METHOD == TEST
   unsigned **ref = create_sand_array(DIM);
   unsigned long ref_time = 0;
   int repeat = 1;

   // NOTE: We use naive compute time for reference
   ref_time = process("SEQ REF",
   		      ref, ref, compute_naive, ref_time, true, repeat);

   ref_time = fmin(ref_time,
   		   process ("SEQ EUCL",
   			    ref, sand, compute_eucl, ref_time, true, repeat));

   ref_time = fmin(ref_time,
   		   process ("SEQ EUCL CHUNK",
   			    ref, sand, compute_eucl_chunk, ref_time,
   			    true, repeat));

   // NOTE: We use best sequential time for reference
   process ("PAR OMP",
   	    ref, sand, compute_omp, ref_time, false, repeat);
   /* process ("PAR OMP TILE", */
   /* 	    ref, sand, compute_omp_tile, ref_time, false, repeat); */
   process ("PAR OMP SEM",
   	    ref, sand, compute_omp_sem, ref_time, false, repeat);

   fprintf(stderr,"\n");

  return 0;
#endif // METHOD TEST
 }
