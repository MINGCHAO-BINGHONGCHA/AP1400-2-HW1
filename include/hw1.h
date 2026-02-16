/*** 
 * @Author: Ming ming_ice_tea@163.com
 * @Date: 2026-02-13 19:06:30
 * @LastEditors: Ming ming_ice_tea@163.com
 * @LastEditTime: 2026-02-16 19:14:16
 * @FilePath: /AP1400-2-HW1/include/hw1.h
 * @Description: 
 * @
 * @Copyright (c) 2026 by Ming ming_ice_tea@163.com, All Rights Reserved. 
 */
#ifndef AP_HW1_H
#define AP_HW1_H

#include <vector>


namespace algebra{
    using Matrix = std::vector<std::vector<double>>;

    Matrix zeros(std::size_t , std::size_t );

    Matrix ones(std::size_t , std::size_t );

    Matrix random(std::size_t , std::size_t , double , double );

    void show(const Matrix& );
     
    Matrix multiply(const Matrix&, double);

    Matrix multiply(const Matrix&, const Matrix&);

    Matrix sum(const Matrix& , double );

    Matrix sum(const Matrix& , const Matrix& );

    Matrix transpose(const Matrix&);

    Matrix minor(const Matrix&, std::size_t , std::size_t );

    double determinant(const Matrix& );

    Matrix inverse(const Matrix& );

    Matrix adjugate(const Matrix& );

    Matrix concatenate(const Matrix&, const Matrix&, int );

    Matrix ero_swap(const Matrix&, std::size_t , std::size_t );
    
    Matrix ero_multiply(const Matrix&, std::size_t, double );

    Matrix ero_sum(const Matrix& , std::size_t , double , std::size_t);

    Matrix upper_triangular(const Matrix&);

}
#endif //AP_HW1_H
