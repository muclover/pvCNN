#pragma once

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>

#include <string>

#include <cryptopp/sha.h>
#include <dirent.h>

#define RED "\033[0;32;31m"
#define NONE "\033[m"
#define YELLOW "\033[1;33m"
#define BLUE "\033[0;34m"

#define log_debug(s) printf(BLUE " ****************** %s ****************** , " YELLOW "LINE:%d , " RED "FILE:%s\n" NONE, s, __LINE__, __FILE__)
#define log_clear() printf(NONE)

std::string string_to_hex(const std::string &in);

std::string SHA256(std::string data);

void SHA256_Bytes(std::string data, unsigned char *output);

template <typename FieldT>
void random_generator(std::vector<FieldT> &a, int n)
{
    srand((unsigned)time(0));
    int l = 0;
    int r = 10;
    for (int i = 0; i < n; ++i)
    {
        a[i] = (rand() % (r - l + 1) + l);
    }
}

template <typename FieldT>
void multi2(const std::vector<FieldT> &a, const std::vector<FieldT> &b, std::vector<FieldT> &c, int M, int N, int K)
{
    for (int i = 0; i < M; ++i) // i: a 行
    {
        for (int j = 0; j < K; ++j) // j: b 列
        {
            for (int k = 0; k < N; ++k) // k为 a 的列 ,b 的行
            {
                c[K * i + j] += a[N * i + k] * b[K * k + j];
            }
        }
    }
}

template <typename T>
void conv_mul(std::vector<T> &in_data, int in_w, int in_h,
              std::vector<T> &kernel, int kernel_w, int kernel_h,
              std::vector<T> &out_data, int out_w, int out_h)
{
    double sum = 0.0;
    for (int i = 0; i < out_h; ++i)
    {
        for (int j = 0; j < out_w; ++j)
        {
            sum = 0.0;
            for (int n = 0; n < kernel_h; ++n)
            {
                for (int m = 0; m < kernel_w; ++m)
                {
                    sum += in_data[(i + n) * in_w + j + m] * kernel[n * kernel_w + m];
                }
            }
            out_data[i * out_w + j] += sum;
        }
    }
}
template <typename T>
void print_vector(const std::vector<T> &v)
{
    for (const auto &e : v)
    {
        std::cout << e << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void print_matrix(const std::vector<T> &matrix, int m, int n)
{
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::cout << matrix[i * n + j] << " ";
        }
        std::cout << ".\n";
    }
}
template <typename T>
void out_matrix_to_file(std::ofstream &out, const std::vector<T> &x, int n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::cout << x[i * n + j] << " ";
            out << x[i * n + j] << " ";
        }
        std::cout << "\n";
        out << "\n";
    }
}

template <typename T>
using Matrix = std::vector<std::vector<T>>;

template <typename T>
void clearVector(std::vector<T> &v)
{
    std::vector<T> vTemp;
    vTemp.swap(v);
}

template <typename T>
T max_element(const Matrix<T> &a)
{
    T max = a[0][0];
    for (auto v : a)
    {
        for (auto ele : v)
        {
            if (ele > max)
            {
                max = ele;
            }
        }
    }
    return max;
}

template <typename T>
T min_element(const Matrix<T> &a)
{
    T min = a[0][0];
    for (auto v : a)
    {
        for (auto ele : v)
        {
            if (ele < min)
            {
                min = ele;
            }
        }
    }
    return min;
}

template <typename T>
std::vector<T> operator*(const Matrix<T> &a, std::vector<T> &b)
{
    std::vector<T> res;

    if (a[0].size() != b.size())
    {
        std::cout << "error" << std::endl;
        throw "error, the col of matrix is not euqal to the size of vector b";
    }

    res.resize(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        for (size_t j = 0; j < b.size(); j++)
        {
            res[i] += a[i][j] * b[j];
        }
    }
    return res;
}

template <typename T>
Matrix<T> operator*(const Matrix<T> &a, const Matrix<T> &b)
{
    Matrix<T> res;
    if (a[0].size() != b.size())
    {
        throw "error in Matrix mul";
    }
    res.resize(a.size());
    for (auto &row : res)
    {
        row.resize(b[0].size());
    }

    for (size_t i = 0; i < a.size(); ++i)
    {
        for (size_t j = 0; j < b[0].size(); ++j)
        {
            for (size_t k = 0; k < b.size(); ++k)
            {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return res;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &m)
{
    for (const auto &row : m)
    {
        for (const auto &ele : row)
        {
            os << ele << " ";
        }
        os << "\n";
    }
    return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
    for (const auto &val : v)
    {
        os << val << " ";
    }
    os << "\n";
    return os;
}
void test_read_dir(std::string PATH);