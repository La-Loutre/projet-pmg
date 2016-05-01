__kernel void sandpiles(__global unsigned *read,
			__global unsigned *write,
			__global unsigned *changed
			)
{
  int x_glob = get_global_id(0);
  int y_glob = get_global_id(1);
  int x_group = get_group_id(0);
  int y_group = get_group_id(1);


  int lx = get_local_id(0);
  int ly = get_local_id(1);
  int k = DIM/TILE; //Integer (16)

  /*
   * x,y real value from 2 dim array
   *
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


  if (y_real > 0 &&
      y_real < DIM-1 &&
      x_real != 0 &&
      x_real != DIM-1)
     {
    int val = read[pos];
    change = (val / 4);
    val %= 4;

    val += read[pos - DIM] / 4
      + read[pos + DIM] / 4
      + read[pos -1] / 4
      + read[pos + 1] / 4;

    write[pos] = val;
    if (change)
      *changed = 1;
   }
  else
    {
      write[pos] = 0;
    }



}
