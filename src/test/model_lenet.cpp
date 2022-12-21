#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "libff/algebra/fields/field_utils.hpp"
#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#include "../qmp/matrix_gadget.hpp"
#include "../qmp/qmp.hpp"
#include "../qmp/timer.hpp"
#include "../qmp/utils.hpp"
#include "../qmp/writer.hpp"

using namespace std;

string weight_file = "/root/pvCNN/src/data/lenet-weight-scale-10.csv";
string input_file = "/root/pvCNN/src/data/lenet-input-scale-10.csv";

int number_of_input = 1;
constexpr int number_of_filter = 1;

template <typename T>
void print_kernel(const vector<vector<vector<vector<T>>>> &kernel, int c_out, int c_in, int k)
{
    for (int i = 0; i < c_in; ++i) {
        for (int j = 0; j < c_out; ++j) {
            for (int p = 0; p < k; ++p) {
                for (int q = 0; q < k; ++q) {
                    cout << kernel[i][j][p][q] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}
void print_input(vector<vector<vector<int>>> input, int c_in, int n)
{
    for (int i = 0; i < c_in; ++i) {
        for (int p = 0; p < n; ++p) {
            for (int q = 0; q < n; ++q) {
                cout << input[i][p][q] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
    cout << endl;
}
// void quan_scale(vector<vector<vector<vector<double>>>> input)
vector<vector<int>> quan_scale(vector<vector<double>> input)
{
    // max
    double max = input[0][0];
    double min = input[0][0];
    for (int i = 0; i < input.size(); ++i) {
        for (int j = 0; j < input[0].size(); ++j) {
            if (input[i][j] > max) {
                max = input[i][j];
            }
            if (input[i][j] < min) {
                min = input[i][j];
            }
        }
    }

    auto scale = (max - min) / 255.0;
    auto zero_point = (255.0 - max / scale);
    vector<vector<int>> res;
    res.reserve(input.size() * input.size());
    for (int i = 0; i < input.size(); ++i) {
        for (int j = 0; j < input[0].size(); ++j) {
            res.emplace_back(round(input[i][j] / scale + zero_point));
        }
    }
    // auto biggest = max_element(input.begin(), input.end());
    // auto smallest = min_element(input.begin(), input.end());
    // cout << "------------------" << endl;
    // cout << *biggest << endl;
    // cout << *smallest << endl;
    // zero_point
    return res;
}

vector<vector<vector<int>>> read_input_from_file(const string &file)
{
    // should read c_in =1
    // repeat 3 times
    int c_in = number_of_input;
    int n = 32;
    vector<vector<vector<int>>> input(c_in, vector<vector<int>>(n, vector<int>(n)));
    ifstream ifs;
    ifs.open(file.c_str());

    if (!ifs.is_open()) {
        throw "can't open file";
    }

    string line;
    getline(ifs, line); // 1 line : 3*32*32
    // while (getline(ifs, line))
    // {
    for (int i = 0; i < c_in; ++i) {
        stringstream ss(line);
        for (int p = 0; p < n; ++p) {
            for (int q = 0; q < n; ++q) {
                ss >> input[i][p][q];
                string s;
                getline(ss, s, ',');
            }
        }
    }
    // }
    // print_input(input, c_in, n);
    return input;
}

// c_in * c_out * k * k
vector<vector<vector<vector<vector<int>>>>> read_weight_from_file(const string &file)
{
    ifstream ifs;
    ifs.open(file.c_str());
    if (!ifs.is_open()) {
        throw "can't open file";
    }

    string line;

    int c_in, c_out, k, batch_size;
    c_out = number_of_filter;
    c_in = 1;
    // c_in = 3;
    batch_size = 1;
    k = 5;

    // vector<vector<vector<vector<int>>>> kernel(c_in, vector<vector<vector<int>>>(c_out, vector<vector<int>>(k, vector<int>(k))));
    vector<vector<vector<vector<vector<int>>>>> kernel(3, vector<vector<vector<vector<int>>>>(c_in, vector<vector<vector<int>>>(c_out, vector<vector<int>>(k, vector<int>(k)))));
    // vector<vector<vector<vector<double>>>> kernel(c_out, vector<vector<vector<double>>>(c_in, vector<vector<double>>(k, vector<double>(k))));

    getline(ifs, line);

    stringstream ss(line);

    for (int c = 0; c < 3; ++c) {
        for (int i = 0; i < c_in; ++i) {
            for (int j = 0; j < c_out; ++j) {
                for (int p = 0; p < k; ++p) {
                    for (int q = 0; q < k; ++q) {
                        ss >> kernel[c][i][j][p][q];
                        string s;
                        getline(ss, s, ',');
                    }
                }
            }
        }
    }
    // print_kernel(kernel, c_out, c_in, k);
    // quan_scale(kernel);
    return kernel;
}
// c_in
// only focus one channel
void process_input(const vector<vector<vector<int>>> &input, int c_in, int n, int k, vector<vector<int>> &Mi)
{
    int block = n - k + 1;
    int N = Mi.size();
    cout << "N=" << N << endl;
    int col = 0;
    c_in = number_of_input;

    for (int i = 0; i < c_in; ++i) {
        vector<int> tmp;
        int index = 0;
        while (index < block) {
            for (int p = 0; p < n; ++p) {
                for (int q = 0; q < k; ++q) {
                    tmp.push_back(input[i][index + q][p]);
                }
            }
            index++;
        }
        // print_vector(tmp);
        for (int t = 0; t < tmp.size(); ++t) {
            Mi[t][i] = tmp[t];
        }
    }
}

// c_in : batch size
// only focus c_in=1
// c_out * c_in * k * k
void process_kernel(const vector<vector<vector<vector<int>>>> &kernel, int c_out, int c_channel, int k, int n, vector<vector<int>> &Mk)
{
    c_channel = 1;
    int block = n - k + 1;
    int N = Mk.size();
    cout << "N= " << N << endl;
    int t = 0;
    // vector<vector<int>> tmp;
    int step = k;
    int row = 0;

    // 1 channel
    // c_out * k * k
    // each proof : c_out*k*k conv 1*n*n
    // output:
    for (int i = 0; i < c_channel; ++i) {
        for (int j = 0; j < c_out; ++j) {
            // vector<vector<int>> tmp;
            // tmp.insert(tmp.end(), kernel[i][j].begin(), kernel[i][j].end());
            vector<int> tmp;
            for (int p = 0; p < k; ++p) {
                for (int q = 0; q < k; ++q) {
                    tmp.push_back(kernel[i][j][q][p]);
                }
            }
            // print_vector(tmp);

            for (int a = 0; a < block; ++a) {
                for (int b = 0; b < block; ++b) {
                    for (int col = 0; col < k * k; ++col) {
                        Mk[row][col + b * step + a * k * n] = tmp[col];
                    }
                    row++;
                }
            }
        }
        cout << endl;
    }
    cout << endl;
}

int main()
{
    int number_execute = 9;
    int var_input[10] = {1, 10, 50, 100, 500, 1000, 2000, 3000, 4000, 4480};
    // int var_input[6] = {500, 1000, 2000, 3000, 4000, 4700};
    while (number_execute < 10) {
        number_of_input = var_input[number_execute];
        number_execute++;

        Writer write("Lenet-layer1.txt");
        Timer timer;

        log_debug("LeNet");

        log_debug("Read Kernel");
        auto kernel1 = read_weight_from_file(weight_file);
        log_debug("Read Input");
        auto input1 = read_input_from_file(input_file);

        log_debug("test read input and kernel");

        int number_of_channel = 3;

        int c_in, c_out, n, k, batch_size;
        c_in = 3;
        n = 32;
        k = 5;
        c_out = number_of_filter;
        batch_size = 1;
        int output_size = n - k + 1; // 28
        int num = 1;
        // not batch size
        vector<vector<vector<int>>> output(c_out, vector<vector<int>>(output_size, vector<int>(output_size, 0)));

        int M1 = n * k * output_size; // 32*5*28 = 4480

        // c_out is the number of filter
        int M2 = c_out * output_size * output_size; // 6*28*28 = 4704
        int M = max(M1, M2);
        cout << "dimension: " << M << endl;
        // left matrix
        // every row:  output_size*output_size = one conv result
        // right matrix
        // each col: unfold every col in input

        // output = Mk * Mi
        vector<vector<int>> Mk(M, vector<int>(M));
        vector<vector<int>> Mi(M, vector<int>(M));

        log_debug("process data");
        write.write_string("\n******************** number of input: " + to_string(number_of_input) + " ****************");
        write.write_string("\n******************** number of filter: " + to_string(number_of_filter) + " ****************");

        vector<vector<vector<int>>> three_result(3, vector<vector<int>>(M, vector<int>(M)));

        for (int count = 0; count < number_of_channel; ++count) {
            write.write_string("\n******************** Channel: " + to_string(count) + " ****************\n");
            process_kernel(kernel1[count], c_out, c_in, k, n, Mk);
            process_input(input1, c_in, n, k, Mi);

            timer.start();
            log_debug("matrix multiply");
            auto res = Mk * Mi;
            timer.stop();
            write.keep("*Matrix Multiply", to_string(timer.elapse_time()));
            timer.clear();

            three_result.push_back(res);

            timer.start();
            vector<int> tmp;
            for (int p = 0; p < res.size(); ++p) {
                for (int q = 0; q < batch_size; ++q) {
                    tmp.push_back(res[p][q]);
                }
            }

            for (int index = 0; index < c_out * output_size * output_size; ++index) {
                for (int i = 0; i < c_out; ++i) {
                    for (int j = 0; j < output_size; ++j) {
                        for (int p = 0; p < output_size; ++p) {
                            output[i][j][p] = tmp[index++];
                        }
                    }
                }
            }
            // timer.stop();
            // write.keep("*output transform", to_string(timer.elapse_time()));
            // timer.clear();

            cout << "********* output *********" << endl;

            log_clear();

            vector<int> a = input1[0][0];

            log_debug("QMP start");
            default_r1cs_gg_ppzksnark_pp::init_public_params();
            typedef libff::Fr<default_r1cs_gg_ppzksnark_pp> FieldT;
            int N = M;

            log_debug("QMP");
            libff::enter_block("QMP phrase Start");
            vector<FieldT> A0(N * N, 1);

            vector<FieldT> A1(N * N);
            vector<FieldT> A2(N * N);
            vector<FieldT> A3(N * N);

            ofstream ofs1("/root/pvCNN/output1.txt");
            ofstream ofs2("/root/pvCNN/output2.txt");
            ofstream ofs3("/root/pvCNN/output3.txt");

            for (int i = 0; i < M; ++i) {
                for (int j = 0; j < M; ++j) {
                    A1[i * M + j] = Mk[i][j];
                    A2[i * M + j] = Mi[i][j];
                    A3[i * M + j] = res[i][j];
                    ofs1 << A1[i * M + j] << " ";
                    ofs2 << A2[i * M + j] << " ";
                    ofs3 << A3[i * M + j] << " ";
                }
                ofs1 << "\n";
                ofs2 << "\n";
                ofs3 << "\n";
            }

            // multi2(A1, A2, A3, N, N, N);
            log_debug("QMP start");

            qmp_zksnark<default_r1cs_gg_ppzksnark_pp> cs = qmp_zksnark<
                default_r1cs_gg_ppzksnark_pp>(move(A0),
                                              move(A1),
                                              move(A2),
                                              move(A3));

            log_debug("CRS Generate");

            libff::enter_block("CRS Phrase");

            timer.start();
            qmp_crs<default_r1cs_gg_ppzksnark_pp> crs =
                qmp_zksnark_generator<default_r1cs_gg_ppzksnark_pp>(cs, N);
            timer.stop();
            write.keep("*QMP CRS", to_string(timer.elapse_time()));
            timer.clear();

            libff::leave_block("CRS Phrase");

            log_debug("Pi Generate");
            timer.start();
            libff::enter_block("Proof Phrase");
            qmp_zksnark_proof<default_r1cs_gg_ppzksnark_pp> pi =
                qmp_zksnark_prover<default_r1cs_gg_ppzksnark_pp>(cs, crs);
            timer.stop();
            write.keep("*QMP Prover", to_string(timer.elapse_time()));
            timer.clear();
            libff::leave_block("Proof Phrase");

            libff::G1_vector<default_r1cs_gg_ppzksnark_pp> g_A_vec(cs.A1.size());
            libff::G2_vector<default_r1cs_gg_ppzksnark_pp> g_B_vec(cs.A2.size());
            libff::G1_vector<default_r1cs_gg_ppzksnark_pp> g_C_vec(cs.A3.size());

            log_debug("auxiliary A1 generate");
            timer.start();
            for (int i = 0; i < cs.A1.size(); ++i) {
                g_A_vec[i] = cs.A1[i] * crs.g;
            }

            for (int i = 0; i < cs.A2.size(); ++i) {
                g_B_vec[i] = cs.A2[i] * crs.h;
            }
            log_debug("auxiliary A2 generate");

            for (int i = 0; i < cs.A3.size(); ++i) {
                g_C_vec[i] = cs.A3[i] * crs.g;
            }
            log_debug("auxiliary A3 generate");
            log_debug("auxiliary information generate");

            auxiliary_input_information<default_r1cs_gg_ppzksnark_pp> aux_input =
                auxiliary_input_information<default_r1cs_gg_ppzksnark_pp>(
                    move(g_A_vec),
                    move(g_B_vec),
                    move(g_C_vec));
            timer.stop();
            write.keep("*QMP Auxiliary Input", to_string(timer.elapse_time()));
            timer.clear();

            log_debug(" Proof Verification ");

            libff::enter_block("Verification phrase");
            timer.start();
            bool result = qmp_zksnark_verifier<default_r1cs_gg_ppzksnark_pp>(crs, pi, aux_input);
            timer.stop();
            write.keep("*QMP Verifier", to_string(timer.elapse_time()));
            timer.clear();
            libff::leave_block("Verification phrase");

            cout << "Verified Status: " << result << endl;

            libff::leave_block("QMP phrase");

            write.write(N);
            write.write_string("-------------------- memory --------------");
            write.write_bit("*QMP CRS size in bits\t", crs.size_in_bits());
            write.write_bit("*QMP Aux size in bits\t", aux_input.size_in_bits());
            write.write_bit("*QMP Proof size in bits\t", pi.size_in_bits());
        }

        timer.start();
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                three_result[0][i][j] = three_result[0][i][j] + three_result[1][i][j] + three_result[2][i][j];
            }
        }
        timer.stop();
        write.write_string("three result add timer: " + to_string(timer.elapse_time()) + " *************");
        timer.clear();
    }
    return 0;
}
