__kernel void sandpiles(__global unsigned *read,
			__global unsigned *write)
{
  int x = get_global_id(0);
  int y = get_global_id(1);
  int lx = get_local_id(0);
  int ly  = get_local_id(1);
  int X = get_group_id(0);
  int Y = get_group_id(1);

  if (y != 0 && y != DIM-1 && x != 0 && x != DIM-1) {
    int val = read[y*DIM+x];
    val &= 3;
    val += read[(y-1)*DIM+x] / 4
      + read[(y+1)*DIM+x] / 4
      + read[y*DIM+x-1] / 4
      + read[y*DIM+x+1] / 4;

    write[y*DIM+x] = val;
  }
}
