#include<stdio.h>
#include "xcl2.hpp"
#include "xclhal2.h"
#include "unistd.h"
#include <vector>
#include<math.h>
#include <errno.h>
#include <thread>

#include <fstream>
#include <string>
#include <chrono>
#include <cmath>
#include "sdx_ker_thread.h" 
#define PRINT_SAMPLE_OUT false
#define ONE_GIG (1024UL*1024UL*1024UL)
using namespace std;


void sdx_ker_thread(cl::Buffer buffer_in1, cl::Buffer buffer_w, cl::Context &context, cl::CommandQueue &q, cl::Kernel &krnl_sdx_cppKernel_top ) {


    cl_int err;


    OCL_CHECK(err, err = krnl_sdx_cppKernel_top.setArg(0,buffer_in1));
    OCL_CHECK(err, err = krnl_sdx_cppKernel_top.setArg(1,buffer_w));
    OCL_CHECK(err, err = krnl_sdx_cppKernel_top.setArg(2,(unsigned int)(NUMBER_OF_DATA_SETS/NUM_CU)));

    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_in1}, 0 /* 0 means from host*/));

    OCL_CHECK(err, err = q.enqueueTask(krnl_sdx_cppKernel_top));

    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_w}, CL_MIGRATE_MEM_OBJECT_HOST));

    OCL_CHECK(err, err = q.finish());
}
