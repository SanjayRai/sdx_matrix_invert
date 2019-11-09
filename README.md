# sdx_matrix_invert
Matrix Invert Kernels for Sdx Platforms

Examples were only built and tested on U250 although most will fit in U200 as well.
For the U200 you can set he NUMBER_OF_KERNELS to a smaller value.
Only Supported Values for NUMBER_OF_KERNELS are 1,2,4.

System requirements:

1. Alveo card (Tested on U200 and U250)
2. Xilinx xVitis 2019.2 tools
3. Xilinx XRT xrt_201920.2.3.1301_7.4.1708-xrt.rpm
4. Xilinx deployment Shell xilinx-u250-xdma-201830.2-2580015.x86_64.rpm (for U250)
4. Xilinx development Shell xilinx-u250-xdma-dev-201830.2-2580015.x86_64.rpm (for U250)

Below is the syntax to build Hardware files on U250

% make all TARGET=hw DEVICE=xilinx_u250_xdma_201830_2

Once built, to Execute the kernel on the U250 Card:

% ./sdx_cppKernel_top ./xclbin/sdx_cppKernel_top.hw.xilinx_u250_xdma_201830_2.xclbin
