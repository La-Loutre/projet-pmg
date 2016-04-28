__kernel void sandpiles(__global unsigned *sand,
			__global unsigned *copy)
{
  int x = get_global_id(0);
  int y = get_global_id(1);
  int lx = get_local_id(0);
  int ly  = get_local_id(1);
  int X = get_group_id(0);
  int Y = get_group_id(1);

  copy[x] = sand[x];
}
