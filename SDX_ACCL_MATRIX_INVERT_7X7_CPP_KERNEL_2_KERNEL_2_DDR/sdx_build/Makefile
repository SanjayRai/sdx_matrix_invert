.PHONY: help

help::
	$(ECHO) "Makefile Usage:"
	$(ECHO) "  make all TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform>"
	$(ECHO) "      Command to generate the design for specified Target and Device."
	$(ECHO) ""
	$(ECHO) "  make clean "
	$(ECHO) "      Command to remove the generated non-hardware files."
	$(ECHO) ""
	$(ECHO) "  make cleanall"
	$(ECHO) "      Command to remove all the generated files."
	$(ECHO) ""
	$(ECHO) "  make build TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform>"
	$(ECHO) "      Command to build xclbin application."
	$(ECHO) ""


DESIGN := sdx_cppKernel_top
TARGETS := hw
TARGET := $(TARGETS)
DEVICE := $(DEVICES)
XCLBIN := ./xclbin

include ./utils.mk

DSA := $(call device2dsa, $(DEVICE))
BUILD_DIR := ./_x.$(TARGET).$(DSA).$(DESIGN)

CXX := $(XILINX_SDX)/bin/xcpp
XOCC := $(XILINX_SDX)/bin/xocc

#Include Libraries
include ../../common_src/opencl/opencl.mk
HOST_SRCS += ../../common_src/xcl2/xcl2.cpp
HOST_HDRS += ../../common_src/xcl2/xcl2.hpp

#CXXFLAGS += -lz -fopenmp
CXXFLAGS += $(opencl_CXXFLAGS) -Wall -O0 -g -std=c++14
CXXFLAGS += -I../../matrix_op_common_src/ -I../src/ -I../../common_src -I../../common_src/xcl2

HOST_SRCS +=  ../../matrix_op_common_src/matrix_operation_wrapper.cpp ../src/tb_matrix_inv.cpp
HOST_HDRS += ../../matrix_op_common_src/matrix_determinant_class.h ../../matrix_op_common_src/matrix_operation_wrapper.h ../../matrix_op_common_src/matrix_utility_class.h ../src/sdx_cppKernel_top.h

KRNL_SRCS += ../../matrix_op_common_src/kernel_WRAPPER.cpp ../../matrix_op_common_src/matrix_operation_wrapper.cpp ../../matrix_op_common_src/sdx_cppKernel_top.cpp
KRNL_HDRS += ../../matrix_op_common_src/matrix_determinant_class.h ../../matrix_op_common_src/matrix_operation_wrapper.h ../../matrix_op_common_src/matrix_utility_class.h ../src/sdx_cppKernel_top.h

# Host compiler global settings
CXXFLAGS += -fmessage-length=0
CXXFLAGS += -Wsign-compare -Wno-unknown-pragmas
LDFLAGS += $(opencl_LDFLAGS) -lrt -lstdc++ 

# Kernel compiler global settings
CLFLAGS += -t $(TARGET) --platform $(DEVICE) --save-temps 
#CLFLAGS += --sp sdx_cppKernel_top_1.m_axi_gmem:DDR[0]  --sp sdx_cppKernel_top_2.m_axi_gmem:DDR[1]  --sp sdx_cppKernel_top_3.m_axi_gmem:DDR[2]  --sp sdx_cppKernel_top_4.m_axi_gmem:DDR[3] 
CLFLAGS += --sp sdx_cppKernel_top_1.m_axi_gmem:DDR[1]  --sp sdx_cppKernel_top_2.m_axi_gmem:DDR[2]
CLFLAGS += -I../../matrix_op_common_src/ -I../src/ -I../../common_src -I../../common_src/xcl2 --kernel_frequency 100



EXECUTABLE = $(DESIGN)

EMCONFIG_DIR = $(XCLBIN)/$(DSA)

BINARY_CONTAINERS += $(XCLBIN)/$(DESIGN).$(TARGET).$(DSA).xclbin
BINARY_CONTAINER_OBJS += $(XCLBIN)/$(DESIGN).$(TARGET).$(DSA).xo

CP = cp -rf

.PHONY: all clean cleanall docs emconfig
all: $(EXECUTABLE) $(BINARY_CONTAINERS) emconfig

.PHONY: exe
exe: $(EXECUTABLE)

.PHONY: build
build: $(BINARY_CONTAINERS)

# Building kernel
$(XCLBIN)/$(DESIGN).$(TARGET).$(DSA).xo:  $(KRNL_SRCS) $(KRNL_HDRS)
	mkdir -p $(XCLBIN)
	$(XOCC) $(CLFLAGS) --temp_dir $(BUILD_DIR) -c -k $(DESIGN) -I'$(<D)' -o'$@' --input_files $(KRNL_SRCS) 
$(XCLBIN)/$(DESIGN).$(TARGET).$(DSA).xclbin: $(BINARY_CONTAINER_OBJS)
	mkdir -p $(XCLBIN)
	$(XOCC) $(CLFLAGS) --temp_dir $(BUILD_DIR) -l $(LDCLFLAGS) --nk $(DESIGN):2 -o'$@' $(+)

# Building Host
$(EXECUTABLE): $(HOST_SRCS) $(HOST_HDRS)
	$(CXX) $(CXXFLAGS) $(HOST_SRCS) -o '$@' $(LDFLAGS)

emconfig:$(EMCONFIG_DIR)/emconfig.json
$(EMCONFIG_DIR)/emconfig.json:
	emconfigutil --platform $(DEVICE) --od $(EMCONFIG_DIR)

check: all
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	$(CP) $(EMCONFIG_DIR)/emconfig.json .
	XCL_EMULATION_MODE=$(TARGET) ./$(EXECUTABLE) $(XCLBIN)/$(DESIGN).$(TARGET).$(DSA).xclbin
else
	 ./$(EXECUTABLE) $(XCLBIN)/$(DESIGN).$(TARGET).$(DSA).xclbin
endif
ifneq ($(TARGET),$(findstring $(TARGET), sw_emu hw hw_emu))
$(warning WARNING:Application supports only sw_emu hw TARGET. Please use the target for running the application)
endif

	sdx_analyze profile -i profile_summary.csv -f html

# Cleaning stuff
clean:
	-rm -rf $(EXECUTABLE) $(XCLBIN)/{*sw_emu*,*hw_emu*} 
	-rm -rf profile_* TempConfig system_estimate.xtxt *.rpt *.csv 
	rm -rf src/*.ll _xocc_* .Xil emconfig.json dltmp* xmltmp* *.log *.jou *.wcfg *.wdb

cleanall: clean
	-rm -rf $(XCLBIN)
	-rm -rf _x.*

