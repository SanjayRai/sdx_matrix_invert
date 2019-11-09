# sdx_matrix_invert
Matrix Invert Kernels for Sdx Platforms

Examples were only built and tested on U250 although most will fit in U200 as well.
For the U200 you can set he NUMBER_OF_KERNELS to a smaller value.
Only Supported Values for NUMBER_OF_KERNELS are 1,2,4.

Below is the syntax to build Hardware files on U250

% make all TARGET=hw DEVICE=xilinx_u250_xdma_201830_2

Once built, to Execute the kernel on the U250 Card:

% ./sdx_cppKernel_top ./xclbin/sdx_cppKernel_top.hw.xilinx_u250_xdma_201830_2.xclbin
