#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#ifdef __APPLE__
#  include <OpenCL/opencl.h>
#else
#  include <CL/opencl.h>
#endif
#include <assert.h>
#include "util.h"
#include "vector.h"
#include "sand.h"

#define MAX_PLATFORMS 3
#define MAX_DEVICES   5

// Application-specific data
#define KERNEL_NAME  "sandpiles"
#define KERNEL_FILE  "src/compute.cl"

unsigned SIZE = DIM;
unsigned TILE = 16;

unsigned *input_data, *output_data,*changed_data,*changed_data_reset;
sand_t ref, sand;

cl_mem input_buffer;  // device memory used for input data
cl_mem output_buffer;  // device memory used for output data
cl_mem changed_buffer;
cl_mem changed_reset_buffer;

static void alloc_buffers_and_user_data(cl_context context)
{
  // CPU side
  input_data = &sand[0][0];
  output_data = malloc(SIZE * SIZE * sizeof(unsigned));
  changed_data = malloc(sizeof(unsigned));
  *changed_data = 0;
  changed_data_reset = malloc(sizeof(unsigned));
  *changed_data_reset = 0;
  // Allocate buffers inside device memory
  input_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(unsigned) * SIZE * SIZE, NULL, NULL);
  if (!input_buffer)
    error("Failed to allocate input buffer");

  output_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(unsigned) * SIZE * SIZE, NULL, NULL);
  if (!output_buffer)
    error("Failed to allocate output buffer");

  changed_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(unsigned) , NULL, NULL);
  if (!changed_buffer)
    error("Failed to allocate changed buffer");

  changed_reset_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(unsigned) , NULL, NULL);
  if (!changed_buffer)
    error("Failed to allocate changed buffer");
}

static void check_output_data(char *name,
			      unsigned long ref_time,
			      unsigned long compute_time)
{
  // XXX:

  int check_matrix(unsigned *reff, unsigned *sand)
  {
    int cpt = 0;
    for(int i = 1; i < DIM-1; i++) {
      for(int j = 1; j < DIM-1; j++) {
	if (reff[i*DIM+j] != sand[i*DIM+j]){
	     cpt++;
	}
      }

    }
    return cpt;
}
  timeandcheck("OPENCL", ref_time, compute_time, &ref[0][0], output_data,
	       (check_func_type)check_matrix);
}

static void free_buffers_and_user_data(void)
{
  free(output_data);
  free(changed_data);
  clReleaseMemObject(input_buffer);
  clReleaseMemObject(output_buffer);
  clReleaseMemObject(changed_buffer);
  clReleaseMemObject(changed_reset_buffer);
}

static void send_input(cl_command_queue queue)
{
  cl_int err;

  err = clEnqueueWriteBuffer(queue, input_buffer, CL_TRUE, 0,
			     sizeof(unsigned) * SIZE * SIZE, input_data, 0, NULL, NULL);
  check(err, "Failed to write to input array");


}

static void retrieve_output(cl_command_queue queue, int mod)
{
  cl_int err;
  if (mod%2 !=0)
    err = clEnqueueReadBuffer(queue, output_buffer, CL_TRUE, 0,
			      sizeof(unsigned) * SIZE * SIZE, output_data, 0, NULL, NULL );
  else
    err = clEnqueueReadBuffer(queue, input_buffer, CL_TRUE, 0,
			      sizeof(unsigned) * SIZE * SIZE, output_data, 0, NULL, NULL );
  check(err, "Failed to read output array");
}
static int is_done(cl_command_queue queue)
{
  cl_int err;

  err = clEnqueueReadBuffer(queue, changed_buffer, CL_TRUE, 0,
			    sizeof(unsigned) ,changed_data, 0, NULL, NULL );

  check(err, "Failed to read changed value");

  //  fprintf(stderr,"CG : %d\n",*changed_data);
  return *changed_data == 0;

}
static void send_reset_changed(cl_command_queue queue)
{
  cl_int err;
  unsigned reset = 0;
  err = clEnqueueWriteBuffer(queue, changed_buffer, CL_TRUE, 0,
			     sizeof(unsigned) , changed_data_reset, 0, NULL, NULL);
  check(err, "Failed to write to changed_reset");

}
void start(sand_t inref, sand_t insand, unsigned long ref_time, bool cpu, bool gpu)
{
  cl_platform_id pf[MAX_PLATFORMS];
  cl_uint nb_platforms = 0;
  cl_int err;                         // error code returned from api calls
  cl_device_type device_type = CL_DEVICE_TYPE_ALL;

  if (gpu && cpu)
    device_type = CL_DEVICE_TYPE_ALL;
  else if (gpu)
    device_type = CL_DEVICE_TYPE_GPU;
  else
    device_type = CL_DEVICE_TYPE_CPU;

  ref = inref;
  sand = insand;

  // Get list of OpenCL platforms detected
  err = clGetPlatformIDs(3, pf, &nb_platforms);
  check(err, "Failed to get platform IDs");


  // For each platform do
  for (cl_int p = 0; p < nb_platforms; p++) {
    cl_uint num;
    int platform_valid = 1;
    char name[1024], vendor[1024];
    cl_device_id devices[MAX_DEVICES];
    cl_uint nb_devices = 0;
    cl_context context;                 // compute context
    cl_program program;                 // compute program
    cl_kernel kernel;

    err = clGetPlatformInfo(pf[p], CL_PLATFORM_NAME, 1024, name, NULL);
    check(err, "Failed to get Platform Info");

    err = clGetPlatformInfo(pf[p], CL_PLATFORM_VENDOR, 1024, vendor, NULL);
    check(err, "Failed to get Platform Info");


    // Get list of devices
    err = clGetDeviceIDs(pf[p], device_type, MAX_DEVICES, devices, &nb_devices);


    if(nb_devices == 0)
      continue;

    // Create compute context with "device_type" devices
    context = clCreateContext (0, nb_devices, devices, NULL, NULL, &err);
    check(err, "Failed to create compute context");

    // Load program source into memory
    const char	*opencl_prog;
    opencl_prog = file_load(KERNEL_FILE);

    // Attach program source to context
    program = clCreateProgramWithSource(context, 1, &opencl_prog, NULL, &err);
    check(err, "Failed to create program");

    // Compile program
    {
      char flags[1024];

      sprintf (flags,
	       "-cl-mad-enable -cl-fast-relaxed-math -DDIM=%d -DSIZE=%d -DTILE=%d -DTYPE=%s",
	       DIM, SIZE, TILE, "unsigned");

      err = clBuildProgram (program, 0, NULL, flags, NULL, NULL);
      if(err != CL_SUCCESS) {
	size_t len;

	// Display compiler log
	clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
	{
	  char buffer[len+1];

	  fprintf(stderr, "--- Compiler log ---\n");
	  clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL);
	  fprintf(stderr, "%s\n", buffer);
	  fprintf(stderr, "--------------------\n");
	}
	if(err != CL_SUCCESS)
	  error("Failed to build program!\n");
      }
    }

    // Create the compute kernel in the program we wish to run
    kernel = clCreateKernel(program, KERNEL_NAME, &err);
    check(err, "Failed to create compute kernel");

    // Allocate and initialize input data
    alloc_buffers_and_user_data(context);

    // Iterate over devices
    for(cl_int dev = 0; dev < nb_devices; dev++) {
      cl_command_queue queue;

      char name[1024];
      cl_device_type dtype;

      err = clGetDeviceInfo(devices[dev], CL_DEVICE_NAME, 1024, name, NULL);
      check(err, "Cannot get type of device");
      err = clGetDeviceInfo(devices[dev], CL_DEVICE_TYPE, sizeof(cl_device_type), &dtype, NULL);
      check(err, "Cannot get type of device");

      // Create a command queue
      queue = clCreateCommandQueue(context, devices[dev], CL_QUEUE_PROFILING_ENABLE, &err);
      check(err,"Failed to create command queue");

      // Write our data set into device buffer
      send_input(queue);
      int swap=0;
      struct timeval t1,t2;
      	double time;
      // Execute kernel
      {
	cl_event prof_event;
	cl_ulong start, end;


	// global domain size for our calculation
	size_t global[2] = { SIZE, SIZE};
	// local domain size for our calculation
	size_t local[2]  = { TILE , TILE};


	// Set kernel arguments
	err = 0;
	err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_buffer);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output_buffer);
	err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &changed_buffer);
	check(err, "Failed to set kernel arguments");



	gettimeofday (&t1, NULL);
	int iter = 1;
	do{
	  err  = clSetKernelArg(kernel, swap, sizeof(cl_mem), &input_buffer);
	  err |= clSetKernelArg(kernel, (swap+1)%2, sizeof(cl_mem), &output_buffer);
	  send_reset_changed(queue);
	err = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global, local,
				     0, NULL, &prof_event);



	check(err, "Failed to execute kernel");
	// Wait for the command commands to get serviced before reading back results
	clFinish(queue);
	swap= (swap+1)%2;

	}while(!(is_done(queue)) );





	clReleaseEvent(prof_event);
      }

      // Read back the results from the device to verify the output
      retrieve_output(queue,swap);

      gettimeofday (&t2,NULL);
      // Check performance
      time= (double)TIME_DIFF(t1, t2);

      // Validate computation
      check_output_data((dtype == CL_DEVICE_TYPE_GPU) ? "OPENCL GPU" : "OPENCL CPU", ref_time,time);



      clReleaseCommandQueue(queue);
    }

    // Cleanup
    free_buffers_and_user_data();

    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseContext(context);
  }
}
