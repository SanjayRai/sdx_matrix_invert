#include "xcl2.hpp"
#include "xclhal2.h"
#include "sdx_cppKernel_top.h" 
#include  "ap_int.h"
#include  "stdint.h"
#ifndef SDX_KER_THREAD_H_
#define SDX_KER_THREAD_H_
void sdx_ker_thread(cl::Buffer buffer_in1, cl::Buffer buffer_w, cl::Context &context, cl::CommandQueue &q, cl::Kernel &krnl_sdx_cppKernel_top);

#endif //SDX_KER_THREAD_H_
