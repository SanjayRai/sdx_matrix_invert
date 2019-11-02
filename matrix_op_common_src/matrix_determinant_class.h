#ifndef MATRIX_DETERMINANT_CLASS_H_
#define MATRIX_DETERMINANT_CLASS_H_

#include <iostream>
#include "ap_fixed.h"

template <class Ti, int DIM>
class matrix_determinant_class { 

public :
    matrix_determinant_class(void){};

    Ti crout_decomposition(Ti A[DIM][DIM]) {
            int i, j, k;
            Ti sum = 1.0e-32;
            Ti det = 1.0;
            Ti L[DIM][DIM];
            Ti U[DIM][DIM];

            for (i = 0; i < DIM; i++) {
                    U[i][i] = 1.0;
            }

            for (j = 0; j < DIM; j++) {
                    for (i = j; i < DIM; i++) {
                            sum = 1.0e-32;
                            for (k = 0; k < j; k++) {
                                    sum = sum + L[i][k] * U[k][j];	
                            }
                            L[i][j] = A[i][j] - sum;
                    }

                    for (i = j; i < DIM; i++) {
                            sum = 1.0e-32;
                            for(k = 0; k < j; k++) {
                                    sum = sum + L[j][k] * U[k][i];
                            }
                            U[j][i] = (A[j][i] - sum) / L[j][j];
                    }
            }
            for (i=0 ; i < DIM; i++) {
                det *= L[i][i];
            }
            return det;
    }

};

#endif //MATRIX_DETERMINANT_CLASS_H_
