__kernel void sandpiles(__global unsigned *read,
			__global unsigned *write,
			__global unsigned *changed
			)
{
  __local int work_tab[TILE+2+2*1][TILE+2+2*1];
  int offset = (2+2*1) /2;
  int x_glob = get_global_id(0);
  int y_glob = get_global_id(1);
  int x_group = get_group_id(0);
  int y_group = get_group_id(1);


  int lx = get_local_id(0);
  int ly = get_local_id(1);


  int k = DIM/TILE; // Integer (16)

  /*
   * x,y real value from 2 dim array
   */
  int x_real = x_group*TILE+lx;
  int y_real = y_group*TILE+ly;
  /*
   * Pos from x,y real value 1 dim array
   */
  int pos = x_real + y_real*k * TILE;

  /*
   * We don't have any change at start
   */
  int change = 0;
  int first_time = 1;
  int iteration=0;
  work_tab[ly+1][lx+1] = read[pos];
  int val = 0;

  if (y_real > 0 &&
      y_real < DIM-1 &&
      x_real != 0 &&
      x_real != DIM-1)
    {

      if (lx == 0)
	work_tab[ly+1][lx] = read[pos -1];
      if (ly == 0)
	work_tab[ly][lx+1] = read[pos - DIM];
      if (lx == TILE-1)
	work_tab[ly+1][lx+2] = read[pos + 1];
      if (ly == TILE-1)
	work_tab[ly+2][lx+1] = read[pos + DIM];
      val = work_tab[ly+1][lx+1];
      barrier(CLK_LOCAL_MEM_FENCE);
      for(int i=0;i<1;++i){

	change = (val / 4);
	val &= 3;

	val +=  work_tab[ly+1+1][lx+1]/ 4
	  + work_tab[ly-1+1][lx+1]/ 4
	  + work_tab[ly+1][lx+1+1]/ 4
	  + work_tab[ly+1][lx-1+1]/ 4;

	if (change)
	  *changed = 1;
      }
      write[pos] = val;

    }
  else {
    write[pos] = 0;
  }
}
