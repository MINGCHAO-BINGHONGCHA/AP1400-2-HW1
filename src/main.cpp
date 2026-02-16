/*** 
 * @Author: Ming ming_ice_tea@163.com
 * @Date: 2026-02-13 19:06:30
 * @LastEditors: Ming ming_ice_tea@163.com
 * @LastEditTime: 2026-02-15 20:46:43
 * @FilePath: /AP1400-2-HW1/src/main.cpp
 * @Description: 
 * @
 * @Copyright (c) 2026 by Ming ming_ice_tea@163.com, All Rights Reserved. 
 */

#include <iostream>
#include <gtest/gtest.h>
#include "hw1.h"

int main(int argc, char **argv)
{
    if (true) // make false to run unit-tests
    {
        // debug section

        algebra::Matrix zero_martix = algebra::zeros(3,3);
        std::cout << "zero_martix" << std::endl;
        algebra::show(zero_martix);

        std::cout << std::endl;

        algebra::Matrix one_martix = algebra::ones(3,3);
        std::cout << "one_martix" << std::endl;
        algebra::show(one_martix);

        std::cout << std::endl;
        algebra::Matrix random_martix = algebra::random(3,3,1.0,20.0);
        std::cout << "random_martix" << std::endl;
        algebra::show(random_martix);

        std::cout << std::endl;
        algebra::Matrix const_multiply_martix = algebra::multiply(random_martix,2);
        std::cout << "const_multiply_martix" << std::endl;
        algebra::show(const_multiply_martix);

        std::cout << std::endl;
        algebra::Matrix martix_multiply_martix = algebra::multiply(random_martix,const_multiply_martix);
        std::cout << "martix_multiply_martix" << std::endl;
        algebra::show(martix_multiply_martix);

        std::cout << std::endl;
        algebra::Matrix const_sum_martix = algebra::sum(random_martix,10.0);
        std::cout << "const_sum_martix" << std::endl;
        algebra::show(const_sum_martix);

        std::cout << std::endl;
        algebra::Matrix martix_sum_martix = algebra::sum(martix_multiply_martix,const_sum_martix);
        std::cout << "martix_sum_martix" << std::endl;
        algebra::show(martix_sum_martix);


        std::cout << std::endl;
        algebra::Matrix tf_martix = algebra::random(2,4,1.0,20.0);
        std::cout << "transpose before martix" << std::endl;
        algebra::show(tf_martix);
        algebra::Matrix transpose_martix = algebra::transpose(tf_martix);
        std::cout << std::endl;
        std::cout << "transpose martix" << std::endl;
        algebra::show(transpose_martix);

        std::cout << std::endl;
        algebra::Matrix tt_martix = algebra::random(3,3,1.0,20.0);
        std::cout << "minor before martix" << std::endl;
        algebra::show(tt_martix);

        std::cout << std::endl;
        algebra::Matrix minor_martix = algebra::minor(tt_martix, 1, 1);
        std::cout << "minor martix" << std::endl;
        algebra::show(minor_martix);

        std::cout << std::endl;
        algebra::Matrix origin_martix = algebra::random(3,3,0.0,10.0);
        std::cout << "orgin_martix" << std::endl;
        algebra::show(origin_martix);
        
        std::cout << "D = " << algebra::determinant(origin_martix) << std::endl;

        std::cout << std::endl;
        
        algebra::Matrix origin_martix1 = {{1.0,2.0,3.0},{0.0,1.0,4.0},{5.0,6.0,0.0}};
        std::cout << "orgin_martix" << std::endl;
        algebra::show(origin_martix1);

        std::cout << std::endl;

        algebra::Matrix inverse_martix = algebra::inverse(origin_martix1);
        std::cout << "inverse martix" << std::endl;
        algebra::show(inverse_martix);
            
        std::cout << std::endl;
        algebra::Matrix origin_matrix2 = algebra::random(3,3,0.0,20.0);
        algebra::Matrix origin_matrix3 = algebra::random(3,3,0.0,20.0);
        std::cout << "orgin_martix1" << std::endl;
        algebra::show(origin_matrix2);
        std::cout << "orgin_martix2" << std::endl;
        algebra::show(origin_matrix3);

        std::cout << std::endl;

        algebra::Matrix concatenate_matrix_zero = algebra::concatenate(origin_matrix2,origin_matrix3,0);
        std::cout << "martix_concatenate_zero" << std::endl;

        algebra::show(concatenate_matrix_zero);

        std::cout << std::endl;

        algebra::Matrix concatenate_matrix_one = algebra::concatenate(origin_matrix2,origin_matrix3,1);
        std::cout << "martix_concatenate_one" << std::endl;

        algebra::show(concatenate_matrix_one);

        std::cout << std::endl;
        std::cout << "orgin_martix" << std::endl;
        algebra::show(origin_martix1);

        algebra::Matrix swap_matrix = algebra::ero_swap(origin_martix1, 0,2);
        std::cout << "swap_matrix" << std::endl;
        algebra::show(swap_matrix);

        std::cout << std::endl;
        std::cout << "orgin_martix" << std::endl;
        algebra::show(origin_martix1);

        algebra::Matrix mu_matrix = algebra::ero_multiply(origin_martix1, 0,2);
        std::cout << "multiply_matrix" << std::endl;
        algebra::show(mu_matrix);

        std::cout << std::endl;
        std::cout << "orgin_martix" << std::endl;
        algebra::show(origin_martix1);

        algebra::Matrix sum_row_matrix = algebra::ero_sum(origin_martix1, 0, 2, 2);
        std::cout << "multiply_matrix" << std::endl;
        algebra::show(sum_row_matrix);

        std::cout << std::endl;
        std::cout << "orgin_martix" << std::endl;
        algebra::show(origin_martix1);

        algebra::Matrix upper_triangular_matrix = algebra::upper_triangular(origin_martix1);
        std::cout << "upper_triangular_matrix" << std::endl;
        algebra::show(upper_triangular_matrix);



    }
    else
    {
        ::testing::InitGoogleTest(&argc, argv);
        std::cout << "RUNNING TESTS ..." << std::endl;
        int ret{RUN_ALL_TESTS()};
        if (!ret)
            std::cout << "<<<SUCCESS>>>" << std::endl;
        else
            std::cout << "FAILED" << std::endl;
    }
    return 0;
}