/*** 
 * @Author: Ming ming_ice_tea@163.com
 * @Date: 2026-02-13 19:06:30
 * @LastEditors: Ming ming_ice_tea@163.com
 * @LastEditTime: 2026-02-17 00:00:52
 * @FilePath: /AP1400-2-HW1/src/hw1.cpp
 * @Description: 
 * @
 * @Copyright (c) 2026 by Ming ming_ice_tea@163.com, All Rights Reserved. 
 */
#include "hw1.h"
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <cstddef>
#include <cmath>

namespace algebra{

    /**
     * @brief 创建一个制定行和列的,所有元素为0的矩阵
     * 
     * @param n 矩阵的行数, n > 0
     * @param m 矩阵的列数, m > 0
     * @return Matrix 返回初始化完成的Matrix;若n,m不合法,返回空的Matrix
     */
    Matrix zeros(size_t n, size_t m){

        if (n < 0 && m < 0){
            return Matrix();
        }

        std::vector<double> row(m,0);
        Matrix matrix(n,row);

        return matrix;

    }

    /**
     * @brief 创建一个制定行和列的,所有元素为1的矩阵
     * 
     * @param n 矩阵的行数, n > 0
     * @param m 矩阵的列数, m > 0
     * @return Matrix 返回初始化完成的Matrix;若n,m不合法,返回空的Matrix
     */
    Matrix ones(size_t n, size_t m){
        
        if (n < 0 && m < 0){
            return Matrix();
        }


        std::vector<double> row(m,1);
        Matrix matrix(n,row);

        return matrix;

    }

    /**
     * @brief 创建一个指定行和列的,元素为指定范围内随机数的矩阵
     * 
     * @param n 矩阵的行数, n >= 0
     * @param m 矩阵的列数, m >= 0
     * @param min 随机数区间的最小值,包含
     * @param max 随机数区间的最大值,不包含
     * @return Matrix 返回初始化完成的Matrix;若n,m不合法,返回空的Matrix;若min > max抛出异常
     */
    Matrix random(size_t n, size_t m, double min, double max){
        
        if(n < 0 && m < 0){
            return Matrix();
        }

        if(min > max){
            throw std::invalid_argument("min > max");
        }

        
        std::mt19937 generator;
        std::uniform_real_distribution<> distribution(min,max);

        
        Matrix matrix;
        for(Matrix::size_type i = 0; i != n; ++i){
            
            std::vector<double> row;

            for(std::vector<double>::size_type j = 0; j != m; ++j){
                row.push_back(distribution(generator));
            }

            matrix.push_back(row);
        }

        return matrix;

    }


    /**
     * @brief 格式化输出矩阵
     * 
     * @param martix 需要输出的矩阵
     */
    
    void show(const Matrix& martix){

        if(martix.empty()){
            std::cerr << "martix is empty" << std::endl;
            return ;
        }
        
        for(Matrix::const_iterator it = martix.begin(); it != martix.end(); ++it){
            for(std::vector<double>::const_iterator it1 = it->begin(); it1 != it->end(); ++it1){
                //std::cout << std::left << std::setw(6) << std::setprecision(3) << *it1 << " ";
                std::cout << std::left << std::setw(6) << std::fixed << std::setprecision(3) << *it1 << " ";
            }
            std::cout << std::endl;
        }

    }
        
    /**
     * @brief 使矩阵乘以一个常数c
     * 
     * @param martix 需要乘法运算的矩阵,不能为空
     * @param c 用来乘矩阵的常数
     * @return Matrix 返回计算完成的Matrix;若矩阵martix不合法,返回空矩阵
     */
    Matrix multiply(const Matrix& martix , double c){
        
        if(martix.empty()){
            return Matrix();
        }

        Matrix multiply_martix(martix);

        for(Matrix::iterator it = multiply_martix.begin(); it != multiply_martix.end(); ++it){
            for(std::vector<double>::iterator it1 = it->begin(); it1 != it->end(); ++it1){
                (*it1) *= c;
            }
            
        }

        return multiply_martix;

    }

    /**
     * @brief 两矩阵相乘
     * 
     * @param matrix1 矩阵1,不能为空
     * @param matrix2 矩阵2,不能为空
     * @return Matrix 返回计算完成的Matrix;若martix1,matrix2不合法,返回空矩阵;matrix1的m != matrix2的n,抛出异常
     */

    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2){
        
        
        if((matrix1.empty()) || (matrix2.empty())){
            return Matrix();
        }

        if(matrix1[0].size() != matrix2.size()){
            throw std::invalid_argument("矩阵维度不匹配");
        }
        

        Matrix multiply_matrix;
        for(int k = 0; k != matrix1[0].size(); ++k){
            std::vector<double> row;
            for(int j=0; j != matrix2[0].size(); ++j){
                
                double sum=0;
                
                for(int i=0; i != matrix1[0].size(); ++i){
                    sum += matrix1[k][i] * matrix2[i][j];
                }
                    
                row.push_back(sum);
            }
            multiply_matrix.push_back(row);
        }

        return multiply_matrix;
        
    }

    /**
     * @brief 将常数c添加到矩阵的每一个元素中
     * 
     * @param matrix 需要计算的矩阵,不能为空
     * @param c 需要计算的常数
     * @return Matrix 返回计算后的Matrix,若matrix不合法,返回空的矩阵
     */
    
    Matrix sum(const Matrix& matrix ,double c){

        if(matrix.empty()){
            return Matrix();
        }

        Matrix sum_matrix;

        for(Matrix::const_iterator it = matrix.begin(); it != matrix.end(); ++it){
            
            std::vector<double> sum;

            for(std::vector<double>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2){
                sum.push_back(c + *it2);
            }

            sum_matrix.push_back(sum);
        }

        return sum_matrix;

    }

    /**
     * @brief 将两个矩阵相加,每个对应的元素位置相加
     * 
     * @param matrix1 矩阵1,不能为空
     * @param matrix2 矩阵2,不能为空
     * @return Matrix 返回计算后的矩阵,若矩阵1和矩阵2为空,返回空矩阵;两矩阵的行列都不相等抛出异常
     */

    Matrix sum(const Matrix& matrix1, const Matrix& matrix2){
        
        if((matrix1.empty()) && (matrix2.empty())){
            return Matrix();
        }
        
        if((matrix1.size() != matrix2.size()) || (matrix1[0].size() != matrix2[0].size())){
            throw std::invalid_argument("矩阵维度不匹配");
        }

        Matrix sum_matrix;
        
        for(Matrix::size_type i = 0; i != matrix1.size(); ++i){
            
            std::vector<double> sum;

            for(std::vector<double>::size_type j = 0; j != matrix1[0].size(); ++j){
                sum.push_back((matrix1[i][j] + matrix2[i][j]));
            }

            sum_matrix.push_back(sum);

        }

        return sum_matrix;

    }

    /**
     * @brief 生成输入矩阵的转置矩阵
     * 
     * @param matrix 需要进行转置的矩阵,不能为空
     * @return Matrix 返回转置后的矩阵;若矩阵为空,返回空矩阵
     */

    Matrix transpose(const Matrix& matrix){

        if(matrix.empty()){
            return Matrix();
        }

        Matrix transpose;

        for(std::vector<double>::size_type i = 0; i != matrix[0].size(); ++i){
            
            std::vector<double> row;

            for(Matrix::size_type j = 0; j != matrix.size(); ++j){
                row.push_back(matrix[j][i]);
            }

            transpose.push_back(row);
        }

        return transpose;
    }

    /**
     * @brief 生成输入矩阵的余子式
     * 
     * @param matrix 需要生成余子式的矩阵,不能为空
     * @param n 对应元素的行的位置, n > 0
     * @param m 对应元素的列的位置, m > 0
     * @return Matrix 返回对应元素的余子式;若矩阵和n,m不合法,返回空矩阵
     */

    Matrix minor(const Matrix& matrix, size_t n, size_t m){
        
        if(n < 0 && m < 0 && (matrix.empty())){
            return Matrix();
        }

        Matrix minor = matrix;

        for(Matrix::iterator it = minor.begin(); it != minor.end(); ++it){
            it->erase(it->begin()+m);
        }

        minor.erase(minor.begin()+n);

        return minor;

    }

    /**
     * @brief 计算输入矩阵的行列式
     * 
     * @param matrix 计算行列式的矩阵,不能为空
     * @return double 返回行列式的值;若矩阵为空,返回1.0;若矩阵的行和列不相等,抛出异常
     */

    double determinant(const Matrix& matrix){

        if(matrix.empty()){
            return 1.0;
        }

        if (matrix.size() != matrix[0].size()){
            throw std::invalid_argument("矩阵不为方阵");
        }

        if(matrix.size() == 1){
            return matrix[0][0];
        }

        double determinant_num = 0;
            
        for(std::vector<double>::size_type i = 0; i != matrix[0].size(); ++i){
            determinant_num += (matrix[0][i] * pow((-1),0+i) * determinant(minor(matrix,0,i)));
        }

        return determinant_num;

    }

    /**
     * @brief 生成输入矩阵的逆矩阵
     * 
     * A^{-1} = 1/det(A) * A^{*}
     * 
     * @param matrix 需要生成逆矩阵的矩阵,不能为空
     * @return Matrix 返回逆矩阵;若输入矩阵为空,返回空矩阵;若矩阵不为方阵或矩阵为奇异矩阵,抛出异常
     */

    Matrix inverse(const Matrix& matrix){
        
        if(matrix.empty()){
            return Matrix();
        }

        if (matrix.size() != matrix[0].size()){
            throw std::invalid_argument("矩阵不为方阵");
        }


        double det = algebra::determinant(matrix);

        if(det == 0){
            throw std::invalid_argument("矩阵为奇异矩阵");
        }
        
        
        Matrix adjugate_matrix = algebra::adjugate(matrix);

        Matrix inverse_matrix = algebra::multiply(adjugate_matrix,(1.0/det));

        return inverse_matrix;

    }


    /**
     * @brief 生成伴随矩阵(A^{*})
     * 
     * @param matrix 需要生成伴随矩阵的矩阵,不能为空
     * @return Matrix 返回伴随矩阵;若输入矩阵为空,返回空矩阵
     */
    Matrix adjugate(const Matrix& matrix){

        if(matrix.empty()){
            return Matrix();
        }

        Matrix Algebraic_cofactor_matrix;

        for(Matrix::size_type i = 0;i != matrix.size(); ++i){

            std::vector<double> row;

            for(std::vector<double>::size_type j = 0;j != matrix[0].size(); ++j){
                row.push_back((pow(-1,i+j) * algebra::determinant(algebra::minor(matrix,i,j))));
            }

            Algebraic_cofactor_matrix.push_back(row);

        }

        Matrix transpose_Algebraic_cofactor_matrix = algebra::transpose(Algebraic_cofactor_matrix);
        
        return transpose_Algebraic_cofactor_matrix;

    }

    /**
     * @brief 根据axis参数,实现两矩阵的堆叠和拼接
     * 
     * @param matrix1 矩阵1,不能为空
     * @param matrix2 矩阵2,不能为空
     * @param axis 参数;0:矩阵堆叠,1:矩阵拼接
     * @return Matrix 返回处理后的矩阵;若两矩阵不合法,返回空矩阵;参数错误,抛出异常
     */
    
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis=0){
        
        if(matrix1.empty() && matrix2.empty()){
            return Matrix();
        }

        if(axis != 0 && axis != 1){
            throw std::invalid_argument("axis error");
        }

        Matrix concatenate_matrix;
        
        if(axis == 0){

            if(matrix1.empty()) return matrix2;
            if(matrix2.empty()) return matrix1;

            if(matrix1[0].size() != matrix2[0].size()){
                throw std::invalid_argument("列错误");
            }

            concatenate_matrix.insert(concatenate_matrix.begin(),matrix1.begin(),matrix1.end());
            concatenate_matrix.insert(concatenate_matrix.begin() + matrix1.size(), matrix2.begin(),matrix2.end());
        }else if(axis == 1){

            if(matrix1.empty()) return matrix2;
            if(matrix2.empty()) return matrix1;

            if(matrix1.size() != matrix2.size()){
                throw std::invalid_argument("行错误");
            }

            concatenate_matrix.resize(matrix1.size());
            for(Matrix::size_type i = 0; i != matrix1.size(); ++i){
                concatenate_matrix[i].insert(concatenate_matrix[i].begin(), matrix1[i].begin(), matrix1[i].end());
            }

            for(Matrix::size_type i = 0; i != matrix2.size(); ++i){
                concatenate_matrix[i].insert(concatenate_matrix[i].end(), matrix2[i].begin(), matrix2[i].end());
            }
            
        }

        return concatenate_matrix;
    }

    /**
     * @brief 矩阵的r1和r2交换
     * 
     * @param matrix 需要交换行的矩阵
     * @param r1 行号1,r1 < matrix.size()
     * @param r2 行号2,r2 < matrix.size()
     * @return Matrix 返回行交换的矩阵;若两矩阵不合法,返回空矩阵;若行号超过了matrix.size(),抛出异常
     */

    Matrix ero_swap(const Matrix& matrix, std::size_t r1, std::size_t r2){
        
        if(matrix.empty()){
            return Matrix();
        }

        if((r1 < 0 || r1 >= matrix.size()) || (r2 < 0 || r2 >= matrix.size()) ){
            throw std::invalid_argument("行号错误");
        }

        Matrix swap_matrix = matrix;
        
        Matrix::iterator temp = swap_matrix.begin() + r1;

        swap_matrix.insert(swap_matrix.begin() + r2, *temp);
        swap_matrix.insert(swap_matrix.begin() + r1, *(swap_matrix.begin() + r2 +1));

        swap_matrix.erase(swap_matrix.begin() + r1 + 1);
        swap_matrix.erase(swap_matrix.begin() + r2 + 1);

        return swap_matrix;

    }


    /**
     * @brief 将第r行的每个元素乘以常数c
     * 
     * @param matrix 矩阵,不能为空
     * @param r 行号,n >= 0
     * @param c 常数
     * @return Matrix 返回指定行乘以常数后的矩阵;若矩阵为空,返回空矩阵;若行号错误,抛出异常
     */
    Matrix ero_multiply(const Matrix& matrix, size_t r, double c){

        if(matrix.empty()){
            return Matrix();
        }

        if(r >= matrix.size() || r < 0 ){
            throw std::invalid_argument("行号错误");
        }

        Matrix multiply_matrix = matrix;

        Matrix::iterator row = multiply_matrix.begin() + r;

        for(std::vector<double>::iterator it = row->begin(); it != row->end(); ++it){
            *(it) *= c;
        }

        return multiply_matrix;

    }

    /**
     * @brief 将矩阵的第r2行加上r1*c
     * 
     * @param matrix 矩阵,不能为空
     * @param r1 需要乘以常数的行的行号 matrix.size() > r1 >= 0
     * @param c  常数
     * @param r2 需要进行加运算行的行号 matrix.size() > r2 >= 0
     * @return Matrix 返回处理后的矩阵;若矩阵为空,返回空矩阵,若行号错误,抛出异常
     */

    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2){

        if(matrix.empty()){
            return Matrix();
        }

        if((r1 < 0 || r1 >= matrix.size()) || (r2 < 0 || r2 >= matrix.size()) ){
            throw std::invalid_argument("行号错误");
        }

        Matrix sum_matrix = ero_multiply(matrix, r1, c);

        for(std::vector<double>::size_type i = 0; i != matrix[0].size(); ++i){
            sum_matrix[r2][i] += sum_matrix[r1][i];
            sum_matrix[r1][i] /= c;
        }

        return sum_matrix;

    }
    

    /**
     * @brief 生成矩阵的上三角矩阵
     * 
     * @param matrix 矩阵, 矩阵不能为空,必须为方阵
     * @return Matrix 返回生成的上三角的矩阵;若矩阵为空,返回空矩阵;矩阵不为方阵,抛出异常
     */

    Matrix upper_triangular(const Matrix& matrix){
        
        if(matrix.empty()){
            return Matrix();
        }

        if(matrix.size() != matrix[0].size()){
            throw std::invalid_argument("矩阵不为方阵");
        }

        Matrix upper_triangular_matrix = matrix;

        size_t rows = upper_triangular_matrix.size();
        size_t cols = upper_triangular_matrix[0].size();

        for(size_t i = 0; i != std::min(rows,cols); ++i){
            
            size_t pivot_row = i;
            while(pivot_row != rows && std::abs(upper_triangular_matrix[pivot_row][i]) < 1e-9){
                pivot_row++;
            }

            if(pivot_row < rows && pivot_row != i){
                upper_triangular_matrix = ero_swap(upper_triangular_matrix, i, pivot_row);
            }

            if(std::abs(upper_triangular_matrix[i][i]) < 1e-9) continue;

            for(size_t j = i+1; j != rows; ++j){
                
                if(std::abs(upper_triangular_matrix[j][i]) > 1e-9){
                    double c = -upper_triangular_matrix[j][i] / upper_triangular_matrix[i][i];
                    upper_triangular_matrix = ero_sum(upper_triangular_matrix, i, c, j);
                }

            }

        }

        return upper_triangular_matrix;
        
    }

}