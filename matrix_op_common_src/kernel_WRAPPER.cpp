#include<stdio.h>
#include "sdx_cppKernel_top.h"

void kernel_WRAPPER (data_t in_arg0[NUM_ELEMENTS_PER_SDX_DATA_BEAT*NUMBER_OF_SDX_BUS_XFERS_PER_INPUT], data_t out_arg0[NUM_ELEMENTS_PER_SDX_DATA_BEAT*NUMBER_OF_SDX_BUS_XFERS_PER_OUTPUT]) {
#pragma HLS PIPELINE II=HLS_KERNEL_WRAPPER_II
data_t fn_in_arg0[DIM][DIM];
#pragma HLS ARRAY_RESHAPE variable=fn_in_arg0 complete dim=0
data_t fn_out_arg0[DIM][DIM];
#pragma HLS ARRAY_RESHAPE variable=fn_out_arg0 complete dim=0
int index;
    for (int row = 0; row < DIM;row++) {
        for (int col = 0; col < DIM;col++) {
            index = col+(row*DIM);
            fn_in_arg0[row][col] = in_arg0[index];
        }
    }

    matrix_operation_wrapper(fn_in_arg0, fn_out_arg0);

    for (int row = 0; row < DIM;row++) {
        for (int col = 0; col < DIM;col++) {
            index = col+(row*DIM);
            out_arg0[index] = fn_out_arg0[row][col];
        }
    }
}
