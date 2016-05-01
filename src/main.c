#include <assert.h>
#include <math.h>
#include <omp.h>
#include <semaphore.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <tgmath.h>
#include <emmintrin.h>
#include <immintrin.h>

#include "display.h"
#include "sand.h"
#include "vector.h"

#define _XOPEN_SOURCE 600

#define TEST 0
#define SEQEUCL 1
#define PAROMP 2
#define PAROMPSEM 3

#define MAX_HEIGHT 4

int compute_chunk(int nthreads)
{
  if ((DIM-2) % nthreads == 0)
    return (DIM-2)/nthreads;
  else
    return 2;
}

static inline int compute_eucl_swap (sand_t sand)
{
  int change = 0;
  sand_t aux = create_sand_array(DIM);

  // We will read the edges, so they should be set to 0
  memset(*aux, 0, DIM*DIM*sizeof(unsigned));

  sand_t swap[2] = {sand, aux};
  sand_t read_from, write_to;
  read_from = swap[0];
  write_to = swap[1];
  int read = 0;
  int write = 1;

  do {
    change = 0;
    for (int y = 1; y < DIM-1; y++) {
      for (int x = 1; x < DIM-1; x++) {
	int val = read_from[y][x];
	// NOTE: works only if MAX_HEIGHT == 4
	change = change | (val >> 2);
	val &= 3;
	val += read_from[y-1][x] / 4
	  + read_from[y+1][x] / 4
	  + read_from[y][x-1] / 4
	  + read_from[y][x+1] / 4;
	write_to[y][x] = val;
      }
    }
    read = 1 - read;
    write = 1 - write;
    read_from = swap[read];
    write_to = swap[write];
  } while(change);

  free(*aux);
  free(aux);
  return change;
}

static inline int compute_eucl_chunk (sand_t sand)
{
  int change = 0;
  int mod4;
  int div4;
  int chunk_size = 2;
  int nb_chunk = (DIM-2) / chunk_size + ((DIM-2)%chunk_size!=0?1:0);
  int chunk[nb_chunk];
  sand_t aux = create_sand_array(DIM);
  sand_t swap[2] = {sand, aux};
  sand_t read_from, write_to;
  read_from = swap[0];
  write_to = swap[1];
  int read = 0;
  int write = 1;

  for (int i=0;i<nb_chunk;++i)
    chunk[i]=1;
  int start=1;
  do{
    change = 0;
  for (int chunk_iter = 0; chunk_iter < nb_chunk; ++chunk_iter) {
    if (!chunk[chunk_iter])
      continue;
    chunk[chunk_iter] = 0;

    for (int y = chunk_iter*chunk_size+start; y < DIM-1 && y < chunk_iter*chunk_size+start+chunk_size; ++y)
      {
	for (int x = 1; x < DIM-1; ++x){
	  int val = read_from[y][x];
	  int save_val = val;
	  // NOTE: works only if MAX_HEIGHT == 4
	  int localchange = val/4;
	  change = change | localchange;
	  val &= 3;
	  val += read_from[y-1][x] / 4
	    + read_from[y+1][x] / 4
	    + read_from[y][x-1] / 4
	    + read_from[y][x+1] / 4;
	  write_to[y][x] = val;

	  chunk[chunk_iter] |= localchange;
	  if (chunk_iter > 0 )
	    chunk[chunk_iter-1] |= (save_val != val);
	  if (chunk_iter  < nb_chunk-1)
	    chunk[chunk_iter+1] |= (save_val != val);

	}
      }


  }
  read = 1 - read;
  write = 1 - write;
  read_from = swap[read];
  write_to = swap[write];
  }while(change);
  free(*aux);
  free(aux);
  return change;
}

static inline int compute_eucl (sand_t sand)
{
  int change = 0;
  int mod4;
  int div4;

  for (int y = 1; y < DIM-1; ++y) {
    for (int x = 1; x < DIM-1; ++x) {
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
	sand[y][x] &= 3;
	sand[y-1][x] += div4;
	sand[y+1][x] += div4;
	sand[y][x-1] += div4;
	sand[y][x+1] += div4;
      }
#endif
    }
  }
  return change;
}

static inline int compute_eucl_vector (sand_t sand)
{
  //FIXME: Only work with multiple of vector size
  int change = 1;
  int mod4;
  int div4;
  __m128i div4_simd,raw_value_simd,mod_simd,previousline_simd,
    nextline_simd,new_value_simd,simd_div,simd_div_mask,simd_shift,simd_shift2,
    simd_mod_mask;
  int value_mod = 3;
  int value_div = 0x3FFFFFFF;
  int value_div_mask[4] ={value_div,
			  value_div,
			  value_div,
			  value_div};
  int value_mod_mask[4] = {value_mod,
			   value_mod,
			   value_mod,
			   value_mod};
  int shift[4] = {2,0,0,0}; //bit inversion
  int shift2[4] = {32,0,0,0};
  int TESTT[4];
  simd_shift = _mm_loadu_si128((__m128i*)&shift);
  simd_shift2 = _mm_loadu_si128((__m128i*)&shift2);
  simd_div_mask = _mm_loadu_si128((__m128i*)&value_div_mask[0]);
  simd_mod_mask = _mm_loadu_si128((__m128i*)&value_mod_mask[0]);

  while(change){
    change = 0;

    for (int y = 1; y < DIM-1; ++y) {
      for (int x = 1; x < DIM-1; x+=4) {
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
	/* SSE part */
	/*
	 * SSE doesn't provide integer
	 * by integer division.
	 * We have to trick
	 * div4=sand[y][x] >> 2; for 4 elements
	 */
	//Copy 4 element into 128b vector (integer)
	raw_value_simd = _mm_loadu_si128((__m128i*)&sand[y][x]);
	//Shift right 2 bits for each .
	// 1 * 8 = 2*4
	div4_simd = _mm_srl_epi32(raw_value_simd,simd_shift);
	// We have to delete overflow by
	// applying a mask (and 00...0011 000...11 000...11)
	div4_simd = _mm_and_si128(div4_simd,simd_div_mask);

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
	previousline_simd = _mm_loadu_si128((__m128i*)&sand[y-1][x]);
	nextline_simd = _mm_loadu_si128((__m128i*)&sand[y+1][x]);
	previousline_simd = _mm_add_epi32(previousline_simd,
					  div4_simd);
	nextline_simd = _mm_add_epi32(nextline_simd,
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
	_mm_storeu_si128((__m128i*)&sand[y-1][x],previousline_simd);
	//Store back next line
	_mm_storeu_si128((__m128i*)&sand[y+1][x],nextline_simd);

	_mm_storeu_si128((__m128i*)&sand[y][x],new_value_simd);
	_mm_storeu_si128((__m128i*)&TESTT,div4_simd);

	sand[y][x-1] += TESTT[0];
	sand[y][x+4] += TESTT[3];
	sand[y][x+1] += TESTT[0] + TESTT[2];
	sand[y][x] += TESTT[1];
	sand[y][x+2] += TESTT[1] + TESTT[3];
	sand[y][x+3] += TESTT[2];
	change = change | TESTT[0] | TESTT[1] | TESTT[2] | TESTT[3];
	/* //First and last element can be done */
	/* //directly (x-1 -> x+4) */
	/* sand[y][x-1] += sand[y][x]; */
	/* sand[y][x+4] += sand[y][x+3]; */

	/* //For others line we must */
	/* // do a swap before */
	/* x1_tmp = sand[y][x+1]; */

	/* sand[y][x+1] += sand[y][x] + sand[y][x+2]; */
	/* sand[y][x] += x1_tmp; */
	/* x2_tmp = sand[y][x+2]; */
	/* sand[y][x+2] += x1_tmp+sand[y][x+3]; */

	/* sand[y][x+3] += x2_tmp; */

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
    int chunk = compute_chunk(nthreads);
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
	  val &= 3;
	  val += sand[y-1][x] / 4
	    + sand[y+1][x] / 4
	    + sand[y][x-1] / 4
	    + sand[y][x+1] / 4;
	  mysand[y][x] = val;

	}
      } // END PARALLEL FOR

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
  } // END PARALLE
  return change;
}

static int iterations = 10;
static inline int compute_omp_iter_v2 (sand_t sand)
{
  int iter_target = iterations;
  int change = 0;
  int nthreads = omp_get_max_threads();
  int chunk = (DIM-2)/nthreads+ ((DIM-2)%nthreads!=0?1:0);
  #pragma omp parallel
  {

    int myid = omp_get_thread_num();
    int mychange = 0;
    sand_t read_buffer = create_sand_array(DIM);
    sand_t write_buffer = create_sand_array(DIM);
    unsigned **swap[2] = {read_buffer,write_buffer};
    int read = 0;
    int write = 1;
    unsigned **read_from = swap[read];
    unsigned **write_to = swap[write];

    do{
      int ysave = 0;
      int my_iter_target = iter_target;
#pragma omp barrier
#pragma omp single // barrier
      change = 0;


#pragma omp for schedule(static, chunk) reduction(|:change)
      /*
	Each thread calculate its chunk
	for 1 iteration.
	For first iteration we read from sand.
      */
      for (int y = 1; y < DIM-1; y++) {
	if (!ysave)
	    ysave = y;
	for (int x = 1; x < DIM-1; x++) {
	  int  val = sand[y][x];
	  change = change | (val /4 );
	  val %= 4;
	  val += sand[y-1][x] / 4
	    + sand[y+1][x] / 4
	    + sand[y][x-1] / 4
	    + sand[y][x+1] / 4;
	  write_to[y][x] = val;
	}
      }
      //Now we have to calculate more lines
      int nb_iter;
      for (nb_iter =0;nb_iter < iter_target-1;++nb_iter){
      	if (nb_iter == 0)
      	  read_from=sand;
      	else{
      	  read = 1 - read;
      	  write = 1 - write;
      	  read_from=swap[read];
      	  write_to=swap[write];
      	}

      	for (int i = 1 ; i < my_iter_target;++i)
      	    {
      	      int Y = ysave-i;
      	      //ieme upper line
      	      if (Y >= 1)
      		{
		  // for(;Y < ysave;++Y)
		    for (int X = 1; X < DIM-1; X++) {

      		    int  val = read_from[Y][X];
      		    //change = change | (val /4 );
      		    val %= 4;
      		    val += read_from[Y-1][X] / 4
      		      + read_from[Y+1][X] / 4
      		      + read_from[Y][X-1] / 4
      		      + read_from[Y][X+1] / 4;
      		    write_to[Y][X] = val;

      		  }

      		}
      	      Y = ysave+chunk+i;
      	      //ieme bottom line
      	      if (Y < DIM-1)
      		{
      		  for (int X = 1; X < DIM-1; X++) {

      		    int  val = read_from[Y][X];
      		    val %= 4;
      		    val += read_from[Y-1][X] / 4
      		      + read_from[Y+1][X] / 4
      		      + read_from[Y][X-1] / 4
      		      + read_from[Y][X+1] / 4;
      		    write_to[Y][X] = val;

      		  }

      		}
      	    }

      	my_iter_target -=1;
      }

      for (int y = ysave; y < chunk+ysave && y < DIM-1; y++) {
      	for (int x = 1; x < DIM-1; x++) {
	  int val;
	  val = write_to[y][x];
      	  sand[y][x] = val;
      	}
      }
      /* if (nb_iter % 2 ==0) */
      /* 	{ */
      /* 	  read = 1 - read; */
      /* 	  write = 1 - write; */
      /* 	  read_from=swap[read]; */
      /* 	  write_to=swap[write]; */
      /* 	} */
      /* if (myid == 0){ */
      /* 	printf("RESULT\n"); */
      /* 	print_matrix(sand,DIM); */
      /* } */
      /* } */

      /* printf("\n"); */


    }while(change);

  }

}
static inline int compute_omp_iter (sand_t sand)
{
  int iter = 1;
  int sup = iter-1;
  int change = 0;

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int myid = omp_get_thread_num();
    int mysup = sup;
    int mychange = 0;

    int chunk = compute_chunk(nthreads);

    unsigned mysand [DIM][DIM];

    do {

#pragma omp barrier
#pragma omp single // barrier
      change = 0;


      // FRONTIER INF
      int inf_frontier = myid*chunk+1;
      for (int y = inf_frontier-sup; y < inf_frontier && y > 1; y++) {
	assert(false);
	for (int x = 1; x < DIM-1; x++) {
	  int val = sand[y][x];
	  val &= 3;
	  val += sand[y-1][x] / 4
	    + sand[y+1][x] / 4
	    + sand[y][x-1] / 4
	    + sand[y][x+1] / 4;
	  mysand[y][x] = val;
	}
      }

      // ITERATIONS
      for (int it = 0; it < iter; it++) {
#pragma omp for schedule(static, chunk)
	for (int y = 1; y < DIM-1; y++) {
	  for (int x = 1; x < DIM-1; x++) {
	    int val = sand[y][x];
	    // NOTE: works only if MAX_HEIGHT == 4
	    mychange = mychange | (val >> 2);
	    val &= 3;
	    val += sand[y-1][x] / 4
	      + sand[y+1][x] / 4
	      + sand[y][x-1] / 4
	      + sand[y][x+1] / 4;
	    mysand[y][x] = val;
	  }
	} // END PARALLEL FOR

	// FRONTIER SUP
	int sup_frontier = myid*chunk+chunk+1;
	for (int y = sup_frontier; y < sup_frontier+sup && y < DIM-1; y++) {
	  assert(false);
	  for (int x = 1; x < DIM-1; x++) {
	    int val = sand[y][x];
	    val &= 3;
	    val += sand[y-1][x] / 4
	      + sand[y+1][x] / 4
	      + sand[y][x-1] / 4
	      + sand[y][x+1] / 4;
	    mysand[y][x] = val;
	  }
	}
	--mysup;
      } // END ITERATIONS

#pragma omp atomic update
      change |= mychange;
#pragma omp barrier
      // SYNCHRONISATION
      if (change) {
#pragma omp for schedule(static, chunk)
	for (int y = 1+sup; y < DIM-1; y++) {
	  for (int x = 1+sup; x < DIM-1; x++) {
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

  int nthreads = omp_get_max_threads();
  int change = 0;
  int tile = 1000;
  int chunk = (DIM-2)/nthreads;
  int nchunk = (DIM-2)/ chunk + ((DIM-2) % chunk == 0 ? 0 : 1);
  int ntile = (DIM-2)/tile +  ((DIM-2) % tile == 0 ? 0 : 1);

#pragma omp parallel shared(change)
  {
    int myid = omp_get_thread_num();
    unsigned mysand [DIM][DIM]; // NOTE: base pointer offset, should be fast

    do {

#pragma omp barrier
#pragma omp single // barrier
      change = 0;

#pragma omp for schedule(static, 1) reduction(|:change)
      for (int c = 0; c < nchunk; c++) {
	for (int t = 0; t < ntile; t++)
	  for (int y = 0; y+c*chunk < DIM-2 && y < chunk; y++)
	    for (int x = 0;  x+t*tile < DIM-2 && x < tile; x++) {
	      int X = x+t*tile+1;
	      int Y = y+c*chunk+1;
	      int val = sand[Y][X];
	      // NOTE: works only if MAX_HEIGHT == 4
	      change = change | (val >> 2);
	      val %= MAX_HEIGHT;
	      val += sand[Y-1][X] / 4
		+ sand[Y+1][X] / 4
		+ sand[Y][X-1] / 4
		+ sand[Y][X+1] / 4;
	      mysand[Y][X] = val;

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

static inline int compute_omp_swap (sand_t sand)
{
  int change = 0;
  sand_t aux = create_sand_array(DIM);

  // We will read the edges, so they should be set to 0
  memset(*aux, 0, DIM*DIM*sizeof(unsigned));

  sand_t swap[2] = {sand, aux};

#pragma omp parallel shared(change)
  {
    sand_t read_from, write_to;
    read_from = swap[0];
    write_to = swap[1];

    int nthreads = omp_get_max_threads();
    int myid = omp_get_thread_num();
    int chunk = compute_chunk(nthreads);
    int read = 0;
    int write = 1;

    do {
#pragma omp barrier
#pragma omp single // barrier
      change = 0;

#pragma omp for schedule(static, chunk) reduction(|:change)
      for (int y = 1; y < DIM-1; y++) {
	for (int x = 1; x < DIM-1; x++) {
	  int val = read_from[y][x];
	  // NOTE: works only if MAX_HEIGHT == 4
	  change = change | (val >> 2);
	  val &= 3 ;
	  val += read_from[y-1][x] / 4
	    + read_from[y+1][x] / 4
	    + read_from[y][x-1] / 4
	    + read_from[y][x+1] / 4;
	  write_to[y][x] = val;
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
  return change;
}

static inline int compute_omp_swap_tile (sand_t sand)
{
  int change = 0;
  int nthreads = omp_get_max_threads();
  sand_t aux = create_sand_array(DIM);
  int tile = 14;
  int chunk = (DIM-2)/nthreads;
  int nchunk = (DIM-2)/ chunk + ((DIM-2) % chunk == 0 ? 0 : 1);
  int ntile = (DIM-2)/tile +  ((DIM-2) % tile == 0 ? 0 : 1);
  // We will read the edges, so they should be set to 0
  memset(*aux, 0, DIM*DIM*sizeof(unsigned));

  sand_t swap[2] = {sand, aux};


#pragma omp parallel shared(change)
  {
    sand_t read_from, write_to;
    read_from = swap[0];
    write_to = swap[1];
    int read = 0;
    int write = 1;

    do {
#pragma omp barrier
#pragma omp single // barrier
      change = 0;

#pragma omp for schedule(static, 1) reduction(|:change)
      for (int c = 0; c < nchunk; c++) {
	for (int t = 0; t < ntile; t++)
	  for (int y = 0; y+c*chunk < DIM-2 && y < chunk; y++)
	    for (int x = 0;  x+t*tile < DIM-2 && x < tile; x++) {
	      int X = x+t*tile+1;
	      int Y = y+c*chunk+1;
	      int val = read_from[Y][X];

	      change = change | (val >> 2);
	      val &= 3;
	      val += read_from[Y-1][X] / 4
		+ read_from[Y+1][X] / 4
		+ read_from[Y][X-1] / 4
		+ read_from[Y][X+1] / 4;

	      write_to[Y][X] = val;
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
  return change;
}

static inline int compute_omp_swap_nowait (sand_t sand)
{
  // XXX: not finished and incorrect

  int nthreads = omp_get_max_threads();
  sand_t aux = create_sand_array(DIM);

  // We will read the edges, so they should be set to 0
  memset(*aux, 0, DIM*DIM*sizeof(unsigned));

  sand_t swap[2] = {sand, aux};

  sem_t *toto = malloc(sizeof(sem_t));
  sem_init(toto, 0, 0);
  sem_t *locks = malloc(sizeof(sem_t)*(nthreads-1));
  for (int i = 0; i < nthreads-1; i++)
    assert(sem_init(&locks[i], 0, 0) ==0);

  int change = nthreads;
  int last_change = 0;
  bool last_iteration = false;

#pragma omp parallel
  {
    sand_t read_from, write_to;
    read_from = swap[0];
    write_to = swap[1];

    int myid = omp_get_thread_num();
    int chunk = compute_chunk(nthreads);
    int read = 0;
    int write = 1;
    int mychange;

    do {
      mychange = 0;
#pragma omp for schedule(static, chunk)
      for (int y = 1; y < DIM-1; y++) {
	int chunk_number = (y-1) / chunk;
	// if nthreads is not a multiple of DIM
	// NOTE: two incrorrect branch predictions at maximum
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
	  mychange = mychange | (val >> 2);
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

      // prepare for next iteration
      read = 1 - read;
      write = 1 - write;
      read_from = swap[read];
      write_to = swap[write];

      if (!last_iteration) {
	// we continue until our chunk is balanced
	if (mychange != 0)
	  continue;

	else {
	  int tmp;
#pragma omp atomic capture
	  tmp = --change;
	  // we wait for the last thread to finish
	  if (tmp != 0)
	    sem_wait(toto);
	  else {
	    // we make a last iteration to be sure nobody can make changes
	    last_iteration = true;
	    for (int i = 0; i < nthreads-1; i++)
	      sem_post(toto);
	  }
	}
      }

      // a last iteration is made to be sure the changes from the last
      // thread had no impact on the rest of the threads
      else {
#pragma omp atomic update
	last_change |= mychange;
#pragma omp barrier
	if (last_change == 0)
	  break; // END DO WHILE
	else {
	  // someone has still some work to do
#pragma omp pragma single // barrier
	  {
	    change = nthreads;
	    last_change = 0;
	    last_iteration = false;
	  }
	  continue;
	}
      }

    } while(true);
  } // END PARALLEL
  free(*aux);
  free(aux);
  free(locks);
  return change;
}

int main (int argc, char **argv)
{

  // omp_set_nested(1);
  //   omp_set_num_threads(4);


  printf("BINDING %d ", omp_get_proc_bind());
  printf("DIM %d CASE %d\n", DIM, CASE);

  sand_t sand = create_sand_array(DIM);
  sand_init (sand);

#if METHOD == SEQEUCL
  display_init (argc, argv,
		DIM,
		MAX_HEIGHT,
		get,
		iterate,
		compute_naive,
		sand);
  return 0;
#endif // METHOD SEQ REF

#if METHOD == PAROMP
  display_init (argc, argv,
		DIM,
		MAX_HEIGHT,
		get,
		iterate,
		compute_eucl_swap,
		sand);
  return 0;
#endif // METHOD SEQ EUCL SWAP

#if METHOD == PAROMPSEM
  display_init (argc, argv,
		DIM,
		MAX_HEIGHT,
		get,
		iterate,
		compute_omp_swap,
		sand);
  return 0;
#endif // METHOD PAR OMP SWAP

#if METHOD == TEST
  sand_t ref = create_sand_array(DIM);
  unsigned long ref_time = 0;
  int repeat = 1;

  // NOTE: We use the previous best compute time for reference

  ref_time = process("SEQ REF",
  		     ref, ref, compute_naive, ref_time, true, repeat);

  ref_time = fmin(ref_time,
  		  process ("SEQ EUCL",
  			   ref, sand, compute_eucl, ref_time, true, repeat));

  ref_time = fmin(ref_time,
  		  process ("SEQ EUCL SWAP",
  			   ref, sand, compute_eucl_swap, ref_time,
  			   true, repeat));

  /* ref_time = fmin(ref_time, */
  /* 		  process ("SEQ EUCL CHUNK", */
  /* 			   ref, sand, compute_eucl_chunk, ref_time, */
  /* 			   false, repeat)); */

  /* ref_time = fmin(ref_time, */
  /* 		  process ("SEQ EUCL VECTOR", */
  /* 			   ref, sand, compute_eucl_vector, ref_time, */
  /* 			   false, repeat)); */

  //  NOTE: We use best sequential time for reference

  int max = omp_get_max_threads();
  for (int i = 1; i <= max; i++) {
    omp_set_num_threads(i);

    printf("NTHREADS %d\n", omp_get_max_threads());

    process ("PAR OMP",
    	     ref, sand, compute_omp, ref_time, false, repeat);

    /* process ("PAR OMP TILE", */
    /* 	     ref, sand, compute_omp_tile, ref_time, false, repeat); */

    process ("PAR OMP SWAP",
  	     ref, sand, compute_omp_swap, ref_time, false, repeat);

    /* process ("PAR OMP SWAP TILE", */
    /* 	     ref, sand, compute_omp_swap_tile, ref_time, false, repeat); */

    /* process ("PAR OMP SWAP NOWAIT", */
    /* 	     ref, sand, compute_omp_swap_nowait, ref_time, false, repeat); */

    /* iterations = 8; */
    /* process ("PAR OMP ITER", */
    /* 	     ref, sand, compute_omp_iter_v2, ref_time, false, repeat); */

  }

  sand_init(sand);
  // OPENCL GPU
  start(ref, sand, ref_time, false, true);
  puts("\n");
  return 0;
#endif // METHOD TEST
}
