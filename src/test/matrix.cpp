/*
  one kernel and one input to prove the generate correct result
*/
#include <iostream>
#include <vector>

#include "libff/algebra/fields/field_utils.hpp"
#include "libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp"
#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>

#include "../qmp/utils.hpp"
#include "../qmp/qmp.hpp"
#include "../qmp/matrix_gadget.hpp"

using namespace libsnark;
using namespace std;

string input_file = "./src/data/M-input_3000.txt";
string filter_file = "./src/data/M-filter_3000.txt";

int main()
{
    default_r1cs_gg_ppzksnark_pp::init_public_params();
    typedef libff::Fr<default_r1cs_gg_ppzksnark_pp> FieldT;

    log_debug("QMP");

    libff::enter_block("QMP phrase");
    float x;
    int c = 1;

    int dimension = 28;
    int weight = 5;
    int block = dimension - weight + 1; // 24
    vector<float> first_input(c * dimension * dimension);

    ifstream ifs(input_file);

    log_debug("读取输入数据");
    //读取输入数据
    for (int i = 0; i < first_input.size(); ++i)
    {
        ifs >> x;
        first_input[i] = x;
    }
    ifs.close();

    int filter_number = 1;
    vector<float> first_filter(filter_number * weight * weight);
    ifs.open(filter_file);
    for (int i = 0; i < first_filter.size(); ++i)
    {
        ifs >> x;
        first_filter[i] = x;
    }
    ifs.close();

    log_debug("proof started");

    log_debug("计算 input 展开的矩阵");
    int M_row = block * weight * dimension;
    int M_column = M_row;

    Matrix<float> M_input(M_row); // matrix
    int d = 0;
    for (int i = 0; i < M_input.size(); ++i)
    {
        M_input[i].resize(M_row);
        for (int j = 0; j < M_row; ++j)
        {
            M_input[i][j] = 0;
        }
    }
    vector<float> input(dimension * dimension);
    int index = 0;
    while (index < c)
    {
        d = 0;
        for (int i = 0; i < dimension; ++i)
        {
            for (int j = 0; j < dimension; ++j)
            {
                input[i * dimension + j] = first_input[(index * dimension + i) * dimension + j];
            }
        }
        for (int i = 0; i < block; ++i)
        {
            for (int j = 0; j < dimension; ++j)
            {
                for (int k = 0; k < weight; ++k)
                {
                    M_input[d][index] = input[dimension * i + k * dimension + j];
                    d++;
                }
            }
        }
        index++;
    }

    vector<float>().swap(input);

    log_debug("计算 weight 展开的矩阵");

    d = 0;
    Matrix<float> M_weight(M_row);
    for (int i = 0; i < M_weight.size(); ++i)
    {
        M_weight[i].resize(M_row);
        for (int j = 0; j < M_row; ++j)
        {
            M_weight[i][j] = 0;
        }
    }

    vector<float> filter(filter_number * weight * weight);
    index = 0;
    while (index < filter_number)
    {
        d = 0;
        for (int i = 0; i < weight; ++i)
        {
            for (int j = 0; j < weight; ++j)
            {
                filter[i * weight + j] = first_filter[(index * weight + i) * weight + j];
            }
        }
        // log_debug("*************");
        // for(int i=0;i<weight;++i){
        //     for(int j=0;j<weight;++j){

        //     }
        // }
        // log_debug("*************");
        for (int t = 0; t < block; ++t)
        {
            for (int i = 0; i < block; ++i)
            {
                for (int j = 0; j < weight; ++j)
                {
                    for (int k = 0; k < weight; ++k)
                    {
                        M_weight[i + t * block + index * block * block][t * weight * dimension + i * weight + d] = filter[k * weight + j];
                    }
                }
                d = 0;
            }
        }
        index++;
    }
    vector<float>().swap(filter);

    Matrix<int> M_first_input_int(M_row);

    for (int i = 0; i < M_first_input_int.size(); ++i)
    {
        M_first_input_int[i].resize(M_row);
        for (int j = 0; j < M_row; ++j)
        {
            M_first_input_int[i][j] = 0;
        }
    }

    log_debug("计算最大最小值");
    auto bigges = max_element<float>(M_input);
    auto smalles = min_element<float>(M_input);

    float scale_in = (bigges - smalles) / 255.0;
    float zero_point_in = round(255.0 - bigges / scale_in);

    for (int i = 0; i < M_input.size(); ++i)
    {
        for (int j = 0; j < M_input.size(); ++j)
        {
            M_first_input_int[i][j] = (round(M_input[i][j] / scale_in) + zero_point_in);
        }
    }
    log_debug("量化weight矩阵");
    Matrix<int> M_first_kernel_int(M_row);
    for (int i = 0; i < M_first_kernel_int.size(); ++i)
    {
        M_first_kernel_int[i].resize(M_row);
        for (int j = 0; j < M_row; ++j)
        {
            M_first_kernel_int[i][j] = 0;
        }
    }
    bigges = max_element<float>(M_weight);
    smalles = min_element<float>(M_weight);

    float scale_kernel = (bigges - smalles) / 255.0;
    float zero_point_kernel = round(255.0 - bigges / scale_kernel);

    for (int i = 0; i < M_weight.size(); ++i)
    {
        for (int j = 0; j < M_weight.size(); ++j)
        {
            M_first_kernel_int[i][j] = (round(M_weight[i][j]) / scale_kernel) + zero_point_kernel;
        }
    }

    log_debug("矩阵相乘");
    Matrix<float> M_mul = M_weight * M_input;
    Matrix<int> M_first_result_int(M_row);
    for (int i = 0; i < M_first_result_int.size(); ++i)
    {
        M_first_result_int[i].resize(M_row);
        for (int j = 0; j < M_row; ++j)
        {
            M_first_result_int[i][j] = 0;
        }
    }

    bigges = max_element<float>(M_mul);
    smalles = min_element<float>(M_mul);
    float scale_result = (bigges - smalles) / 255;
    float zero_point_result = round(255 - bigges / scale_result);
    for (int i = 0; i < M_mul.size(); i++)
    {
        for (int j = 0; j < M_mul.size(); j++)
        {
            M_first_result_int[i][j] = (round((M_mul[i][j]) / scale_result) + zero_point_result);
        }
    }
    float M_scale = scale_kernel * scale_in / scale_result;
    //量化后求解，产生中间量化矩阵
    Matrix<int> M_tempA(M_row);
    //初始化
    for (int i = 0; i < M_tempA.size(); i++)
    {
        M_tempA[i].resize(M_row);
        for (int j = 0; j < M_row; j++)
        {
            M_tempA[i][j] = 0;
        }
    }
    Matrix<int> M_tempB(M_row);
    //初始化
    for (int i = 0; i < M_tempB.size(); i++)
    {
        M_tempB[i].resize(M_row);
        for (int j = 0; j < M_row; j++)
        {
            M_tempB[i][j] = 0;
        }
    }

    for (int i = 0; i < M_first_input_int.size(); i++)
    {
        for (int j = 0; j < M_first_input_int.size(); j++)
        {
            M_tempA[i][j] = M_first_input_int[i][j] - zero_point_in;
            M_tempB[i][j] = M_first_kernel_int[i][j] - zero_point_kernel;
        }
    }

    cout << "矩阵相乘：" << endl;

    Matrix<int> M_result_temp = M_tempB * M_tempA;

    libff::leave_block("QMP phrase");

    log_clear();
}