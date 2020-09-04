#include <stdio.h>
#include <time.h>
#include "CL/opencl.h"

#define NUM_ENTRIES 1024

int main() // (int argc, const char* argv[])
{
  // CONSTANTS
  // The source code of the kernel is represented as a string
  // located inside file: "fft1D_1024_kernel_src.cl". For the details see the
  // next listing.
  const char* KernelSource =
#include "fft.cl"
      ;

  // Looking up the available GPUs
  const cl_uint num = 1;
  clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 0, NULL, (cl_uint*)&num);

  cl_device_id devices[1];
  clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, num, devices, NULL);

  // create a compute context with GPU device
  cl_context context =
      clCreateContextFromType(NULL, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);

  // create a command queue
  clGetDeviceIDs(NULL, CL_DEVICE_TYPE_DEFAULT, 1, devices, NULL);
  cl_command_queue queue = clCreateCommandQueue(context, devices[0], 0, NULL);

  // allocate the buffer memory objects
  cl_mem memobjs[] = {
      clCreateBuffer(
          context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
          sizeof(float) * 2 * NUM_ENTRIES, NULL, NULL),
      clCreateBuffer(
          context, CL_MEM_READ_WRITE, sizeof(float) * 2 * NUM_ENTRIES, NULL,
          NULL)};
  // cl_mem memobjs[0] = // FIXED, SEE ABOVE
  // cl_mem memobjs[1] = // FIXED, SEE ABOVE

  // create the compute program
  // const char* fft1D_1024_kernel_src[1] = {  };
  cl_program program = clCreateProgramWithSource(
      context, 1, (const char**)&KernelSource, NULL, NULL);

  // build the compute program executable
  clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

  // create the compute kernel
  cl_kernel kernel = clCreateKernel(program, "fft1D_1024", NULL);

  // set the args values

  size_t local_work_size[1] = {256};

  clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&memobjs[0]);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&memobjs[1]);
  clSetKernelArg(
      kernel, 2, sizeof(float) * (local_work_size[0] + 1) * 16, NULL);
  clSetKernelArg(
      kernel, 3, sizeof(float) * (local_work_size[0] + 1) * 16, NULL);

  // create N-D range object with work-item dimensions and execute kernel
  size_t global_work_size[1] = {256};

  global_work_size[0] = NUM_ENTRIES;
  local_work_size[0] = 64; // Nvidia: 192 or 256
  clEnqueueNDRangeKernel(
      queue, kernel, 1, NULL, global_work_size, local_work_size, 0, NULL, NULL);
}
