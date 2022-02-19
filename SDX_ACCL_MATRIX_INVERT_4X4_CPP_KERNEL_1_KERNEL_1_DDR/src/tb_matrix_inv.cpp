#include <stdio.h>
#include "xcl2.hpp"
#include "xclhal2.h"
#include "unistd.h"
#include <vector>
#include <math.h>
#include <errno.h>
#include <thread>

#include <fstream>
#include <string>
#include <chrono>
#include <cmath>
#include "sdx_ker_thread.h" 
#include "sdx_cppKernel_top.h" 
#define PRINT_SAMPLE_OUT false
#define ONE_GIG (1024UL*1024UL*1024UL)
using namespace std;

int check_if_Indentity_Matrix (data_t matrix_A[DIM][DIM]){
int id_error_count;
    id_error_count = 0;
    for (unsigned int row = 0; row < DIM;row++) {
        for (unsigned int col = 0; col < DIM;col++) {
            if (row == col) {
                if (fabs(matrix_A[row][col] -1.0f) > ZERO_f){
                    cout << "ERROR DBG :: " << fabs(matrix_A[row][col] -1.0f) << endl;
                    id_error_count++;
                } 
            } else if (fabs(matrix_A[row][col]) > ZERO_f){
                    id_error_count++;
            }
        }
    }
    return id_error_count;
}

void matrix_mult (data_t matrix_A[DIM][DIM], data_t matrix_B[DIM][DIM], data_t mult[DIM][DIM]) {
    for (unsigned int row = 0; row < DIM;row++) {
        for (unsigned int col = 0; col < DIM;col++) {
            mult[row][col] = 0.0f;
        }
    }
    for (unsigned int row = 0; row < DIM;row++) {
        for (unsigned int col = 0; col < DIM;col++) {
          for (unsigned int k = 0; k < DIM;k++) {
            mult[row][col] += matrix_A[row][k]* matrix_B[k][col];
          }
        }
    }
}

void print_matrix (data_t matrix_val[DIM][DIM]) {
    for (unsigned int row = 0; row < DIM;row++) {
        for (unsigned int col = 0; col < DIM;col++) {
            if (fabs(matrix_val[row][col]) < 1.0e-4) {
                matrix_val[row][col] = 0.0f;
            }
            cout << matrix_val[row][col] << " "; 
        }
        cout << endl;
    }
        cout << "------------------------------------"<< endl;

}

void gen_zero_matrix(srai_mem_conv *a) {

    for (unsigned int j = 0 ; j < NUMBER_OF_DATA_SETS; j++) {
        for (unsigned int i = 0 ; i < SDX_CU_LOCAL_IN_SIZE; i++) {
            for (unsigned int row = 0; row < DIM;row++) {
                for (unsigned int col = 0; col < DIM;col++) {
                    a->my_data_t[row][col] = 0.0f; 
                }
            }
            a++;
        }
    }
}
void gen_test_matrix(srai_mem_conv *a) {
//   static data_t temp[DIM][DIM] = {{1.0f, 2.0f, 3.0f, 4.0f},
//                                   {5.0f, 6.0f, 7.0f, 8.0f},
//                                   {2.0f, 6.0f, 4.0f, 8.0f},
//                                   {3.0f, 1.0f, 1.0f, 2.0f}};
   static data_t temp[DIM][DIM] = {{1.0f, ZERO_f, ZERO_f, ZERO_f},
                                   {ZERO_f, 1.0f, ZERO_f, ZERO_f},
                                   {ZERO_f, ZERO_f, 1.0f, ZERO_f},
                                   {ZERO_f, ZERO_f, ZERO_f, 1.0f}};
//    static data_t temp[DIM][DIM] = {{1.0f, 2.0f, 3.0f, 4.0f},
//                                    {5.0f, 6.0f, 7.0f, 8.0f},
//                                    {2.0f, 6.0f, 4.0f, 8.0f},
//                                    {3.0f, 1.0f, 1.0f, 2.0f}};
 
//   static data_t temp[DIM][DIM] = {{1.0f, 0.0f, 0.0f, 0.0f},
//                                   {0.0f, 1.0f, 0.0f, 0.0f},
//                                   {0.0f, 0.0f, 1.0f, 0.0f},
//                                   {0.0f, 0.0f, 0.0f, 1.0f}};
 
//   static data_t temp[DIM][DIM] = {{1.0f, 2.0f, 3.0f, 4.0f},
//                                   {5.0f, 6.0f, 7.0f, 8.0f},
//                                   {9.0f, 10.0f, 11.0f, 12.0f},
//                                   {13.0f, 14.0f, 15.0f, 16.0f}};

  data_t random_scale;

    for (unsigned int j = 0 ; j < NUMBER_OF_DATA_SETS; j++) {
        for (unsigned int i = 0 ; i < SDX_CU_LOCAL_IN_SIZE; i++) {
            //random_scale = (data_t)(rand() % (int)32768.0)/321.01;
            random_scale = (data_t)(1.0f + ((rand() % (int)3276.0)/321.01));
            for (unsigned int row = 0; row < DIM;row++) {
                for (unsigned int col = 0; col < DIM;col++) {
                    //a->my_data_t[row][col] = (temp[row][col])*random_scale; 
                    a->my_data_t[row][col] = (data_t)((rand() % (int)32768.0)/1311.01);
                }
            }
            a++;
        }
    }
}

void print_gen_test_matrix(srai_mem_conv *a) {

    for (unsigned int j = 0 ; j < NUMBER_OF_DATA_SETS; j++) {
        for (unsigned int i = 0 ; i < SDX_CU_LOCAL_IN_SIZE; i++) {
            for (unsigned int row = 0; row < DIM;row++) {
                for (unsigned int col = 0; col < DIM;col++) {
                    cout << a->my_data_t[row][col] << " ";
                }
                cout << endl;
            }
                cout << endl;
                cout << "------------------------------------"<< endl;
                a++;
        }
    }
}

int main(int argc, char** argv) {

    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <XCLBIN File>" << endl;
        return EXIT_FAILURE;
    }

    string binaryFile = argv[1];
    cl_int err;

    time_t t;
    srand((unsigned) time(&t));
    double high_res_elapsed_time = 0.0f;
    double high_res_elapsed_time_HW = 0.0f;
    double high_res_elapsed_time_SW = 0.0f;
    chrono::high_resolution_clock::time_point start_t;
    chrono::high_resolution_clock::time_point stop_t;
    chrono::duration<double> elapsed_hi_res;

  
    cout << "Srai_ DBG NUMBER_OF_DATA_SETS  =  " << NUMBER_OF_DATA_SETS << endl;
    cout << "Srai_ DBG GLOBAL_DATA_IN_SIZE  =  " << GLOBAL_DATA_IN_SIZE << endl;
    cout << "Srai_ DBG GLOBAL_DATA_IN_SIZE_BYTES =  " << GLOBAL_DATA_OUT_SIZE_BYTES << endl;
    if ((GLOBAL_DATA_IN_SIZE_BYTES > (4*ONE_GIG)) | (GLOBAL_DATA_OUT_SIZE_BYTES > (4*ONE_GIG))) {
      cout << "Memory reguirement over 4GB ..........!! \n";
      exit(0);
    }


  
    vector<vector<srai_mem_conv, aligned_allocator<srai_mem_conv>> > a_in_ptr_cx(NUM_CU, vector<srai_mem_conv, aligned_allocator<srai_mem_conv>>(GLOBAL_DATA_IN_SIZE));
    vector<vector<srai_mem_conv, aligned_allocator<srai_mem_conv>> > y_out_ptr_cx(NUM_CU, vector<srai_mem_conv, aligned_allocator<srai_mem_conv>>(GLOBAL_DATA_IN_SIZE));
    //vector<vector<srai_mem_conv, aligned_allocator<srai_mem_conv>> > a_in_ptr_cx(GLOBAL_DATA_IN_SIZE);
    //vector<vector<srai_mem_conv, aligned_allocator<srai_mem_conv>> > y_out_ptr_cx(GLOBAL_DATA_IN_SIZE);

  
    printf("-------------------------------------------------------------\n");
    printf("Create Test Data Set\n");
    printf("Note DATA_IN_SIZE (Input Memory size in bytes  ) = %d (%x)\n",(GLOBAL_DATA_IN_SIZE_BYTES),(GLOBAL_DATA_IN_SIZE_BYTES));
    printf("Note DATA_OUT_SIZE(Input Memory size in bytes  ) = %d (%x)\n",(GLOBAL_DATA_OUT_SIZE_BYTES),(GLOBAL_DATA_OUT_SIZE_BYTES));
    cout << "Size of data_t = " << sizeof(data_t) <<  " Bytes" << endl;
    cout << "Size of srai_mem_conv = " << sizeof(srai_mem_conv) <<  " Bytes" << endl;
    cout << "True Size (in Bytes) of Input Data  = " << sizeof(srai_mem_conv)*GLOBAL_DATA_IN_SIZE<< endl;
    cout << dec;
    printf("-------------------------------------------------------------\n\n\n");

    //Fill ddr4_Memory wr_data_buffer
    cout << "Initializing Memory with InputA args\n";
    for (int i = 0 ; i < NUM_CU; i++) {
        gen_test_matrix(a_in_ptr_cx[i].data());
        //gen_zero_matrix(y_out_ptr_cx[i].data()); // Initialize results pointer with same data as input
    }

    //print_gen_test_matrix(a_in_ptr_c);
    cout << "Memory Initialized with test Data\n";


    //Create Program and Kernel
    vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[0];

    OCL_CHECK(err, cl::Context context(device, NULL, NULL, NULL, &err));


    auto fileBuf = xcl::read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};

    devices.resize(1);
    OCL_CHECK(err, cl::Program program(context, devices, bins, NULL, &err));

    std::string cu_id;
    std::string krnl_name[NUM_CU];
    vector<cl::Kernel> krnl_sdx_cppKernel_top(NUM_CU);
    vector<cl::CommandQueue> q(NUM_CU);
    vector<cl::Buffer> buffer_in1(NUM_CU);
    vector<cl::Buffer> buffer_w(NUM_CU);

    for (int i = 0 ; i < NUM_CU; i++) {
        cu_id = std::to_string(i+1);
        krnl_name[i] = "sdx_cppKernel_top:{sdx_cppKernel_top_" + cu_id + "}";
        OCL_CHECK(err, krnl_sdx_cppKernel_top[i] = cl::Kernel(program,krnl_name[i].c_str(), &err));
        OCL_CHECK(err, q[i] = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));
        OCL_CHECK(err, buffer_in1[i] = cl::Buffer(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, (CHUNK_SIZE*sizeof(srai_mem_conv)), (a_in_ptr_cx[i].data()), &err));
        OCL_CHECK(err, buffer_w[i]   = cl::Buffer(context,CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, (CHUNK_SIZE*sizeof(srai_mem_conv)), (y_out_ptr_cx[i].data()), &err));

    }
  

    cout << "__SRAI DBG :: GLOBAL_DATA_IN_SIZE_BYTES = " << GLOBAL_DATA_IN_SIZE_BYTES << endl;
    cout << "__SRAI DBG :: GLOBAL_DATA_IN_SIZE= " << GLOBAL_DATA_IN_SIZE << endl;
    cout << "__SRAI DBG :: NUMBER_OF_DATA_SETS= " << NUMBER_OF_DATA_SETS << endl;
    cout << "__SRAI DBG :: NUM_CU= " << NUM_CU << endl;
    cout << "__SRAI DBG :: CHUNK_SIZE= " << CHUNK_SIZE << endl;
    cout << "__SRAI DBG :: Chuck_Sz per Kernel GLOBAL_DATA_IN_SIZE_BYTES = " << (CHUNK_SIZE*sizeof(srai_mem_conv)) << endl;
    cout << "__SRAI DBG :: Total Number of " << DIM << "x" << DIM << " Matrices : " <<  (NUMBER_OF_DATA_SETS*SDX_CU_LOCAL_SIZE) << endl;


    thread sdxThread[NUM_CU];
    start_t = chrono::high_resolution_clock::now();
    for (int i = 0 ; i < NUM_CU; i++) {
        sdxThread[i] = thread(sdx_ker_thread, ref(buffer_in1[i]), ref(buffer_w[i]), ref(context), ref(q[i]), ref(krnl_sdx_cppKernel_top[i]));
    }

    for (int i = 0 ; i < NUM_CU; i++) {
        sdxThread[i].join();
    }

    
    stop_t = chrono::high_resolution_clock::now();
    elapsed_hi_res = stop_t - start_t ;
    high_res_elapsed_time_HW = elapsed_hi_res.count();


    printf ("\n-----------------------------------------------------------------\n");
    printf (" ----- Kernel Execution Done Performing Software validation -------\n");
    printf (" ------------   Results  ------------------------------------------\n");
    printf (" ------------------------------------------------------------------\n");
    data_t temp_A[DIM][DIM];
    data_t temp_Y[DIM][DIM];
    data_t result_sw[DIM][DIM];
    data_t mult[DIM][DIM];
    int total_id_error;
    int curr_test_error;
    total_id_error = 0;
    int total_sucess = 0;
    uint32_t random_index[4];

    float cmp_err;
    uint64_t Non_number_count = 0;
    uint64_t singular_matrix_count = 0;
    bool singular_matrix = false;

    srai_mem_conv *a_in_ptr_c[NUM_CU];
    srai_mem_conv *y_out_ptr_c[NUM_CU];

    for (unsigned int i = 0; i < NUM_CU; i++) {
        a_in_ptr_c[i] = a_in_ptr_cx[i].data();
        y_out_ptr_c[i] = y_out_ptr_cx[i].data();
    }
  

    for (unsigned int i = 0; i < 4; i++) {
        random_index[i]  = (uint32_t)(rand() % NUMBER_OF_DATA_SETS); 
    }
    high_res_elapsed_time = 0.0f;
    for (unsigned int cu = 0 ; cu < NUM_CU; cu++) {
        for (unsigned int j = 0 ; j < NUMBER_OF_DATA_SETS/NUM_CU; j++) {
        curr_test_error = 0;
            for (unsigned int i = 0 ; i < SDX_CU_LOCAL_SIZE; i++) {
                for (unsigned int row = 0; row < DIM;row++) {
                    for (unsigned int col = 0; col < DIM;col++) {
                        temp_A[row][col] = a_in_ptr_c[cu]->my_data_t[row][col];
                        temp_Y[row][col] = y_out_ptr_c[cu]->my_data_t[row][col];
                    }
                }
                a_in_ptr_c[cu]++;
                y_out_ptr_c[cu]++;
                matrix_mult (temp_A, temp_Y, mult);
                start_t = chrono::high_resolution_clock::now();
                matrix_operation_wrapper(temp_A, result_sw);
                stop_t = chrono::high_resolution_clock::now();
                elapsed_hi_res = stop_t - start_t ;
                high_res_elapsed_time += elapsed_hi_res.count();
                singular_matrix = false;
                for (unsigned int row = 0; row < DIM;row++) {
                    for (unsigned int col = 0; col < DIM;col++) {
                        if (isnan(temp_Y[row][col])) {
                            Non_number_count++;
                            singular_matrix = true;
                        } else {
                            cmp_err = (fabs(result_sw[row][col]- temp_Y[row][col]));
                            if (cmp_err > 1.0e-6) {
                                curr_test_error++;
                                //cout << "ERR # " << dec << curr_test_error << " :  delta = " << cmp_err << endl;
                                //cout << "Error !! Expected :" << result_sw[row][col] << "  but got : " << temp_Y[row][col] << "  at row = " << row  << " col = " << col << endl;
                            }
                        }
                    }
                }
                if (singular_matrix) {
                    singular_matrix_count++;
                } else if (curr_test_error != 0) {
                    total_id_error++;
    //                cout << "!!Failed Data set : " << j << endl;
    //                   cout << " -------- Expected-------\n";
    //                   print_matrix(result_sw);
    //                   cout << " -------- Got from HW ---\n";
    //                   print_matrix(temp_Y);
                } else {
                    //cout << "  Passed Data set : " << j << endl;
                    total_sucess++;
                }

                if ( ((j == random_index[0]) | (j == random_index[1]) | (j == random_index[2]) | (j == random_index[3])) & PRINT_SAMPLE_OUT ) {
                    cout << "Input Matrix .........\n";
                    print_matrix(temp_A);
                    cout << "Output  Matrix .........\n";
                    print_matrix(temp_Y);
                    cout << "Expected Output  Matrix .........\n";
                    print_matrix(result_sw);
                    cout << "Identity  Matrix .........\n";
                    print_matrix(mult);
                }
            }
        }
    }
    high_res_elapsed_time_SW = high_res_elapsed_time;
    cout << "\n.....................................................\n";
    cout << "SW Execution time =  " <<  high_res_elapsed_time_SW << "s\n";
    cout << "SW THroughput  =  " <<  (GLOBAL_DATA_OUT_SIZE_BYTES/high_res_elapsed_time_SW) << " Bytes/s\n";
    cout << "HW Execution time =  " <<  high_res_elapsed_time_HW << "s\n";
    cout << "HW THroughput =  " <<  (GLOBAL_DATA_OUT_SIZE_BYTES/high_res_elapsed_time_HW) << " Bytes/s\n";
    cout << "Gain (SW_time/HW_time)  =  " << (double)(high_res_elapsed_time_SW/high_res_elapsed_time_HW) << endl;
    printf ("%lu :: Software Number of %lux%lu floating point Matrix Inverts per sec = %f \n",(NUMBER_OF_DATA_SETS*SDX_CU_LOCAL_SIZE), DIM, DIM, (NUMBER_OF_DATA_SETS*SDX_CU_LOCAL_SIZE)/high_res_elapsed_time_SW);
    printf ("%lu :: Hawdware Number of %lux%lu floating point Matrix Inverts per sec = %f \n",(NUMBER_OF_DATA_SETS*SDX_CU_LOCAL_SIZE), DIM, DIM, (NUMBER_OF_DATA_SETS*SDX_CU_LOCAL_SIZE)/high_res_elapsed_time_HW);

    if (total_id_error) {
        cout << "!!!!!!!!!!!!!  TEST UnSucessful! : Total Errors  =  " << total_id_error << endl; 
    } else {
        cout << "TEST Sucessful : Total Errors  =  " << total_id_error << endl; 
    }


    printf (" ------------------------------------------------------------------\n");

    cout << "Results verifcation complete " << endl;
    cout << "Total Error = " << total_id_error << endl;
    cout << "Total Sucess = " << total_sucess << endl;
    cout << "Singular Martix count   = " << singular_matrix_count <<" Percentage of Singular Matrices = " << (100.0*singular_matrix_count)/(NUMBER_OF_DATA_SETS*SDX_CU_LOCAL_SIZE) << "%"<< endl;
    cout << "Total Non numbers (NAN) = " << Non_number_count << endl;

    printf (" ------------------------------------------------------------------\n");



    return 0;
}
