/*
only for matrix multiplication use QMP
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <unistd.h>

#include "libff/algebra/fields/field_utils.hpp"
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>

#include "../qmp/utils.hpp"
#include "../qmp/matrix_gadget.hpp"
#include "../qmp/qmp.hpp"
#include "../qmp/timer.hpp"
#include "../qmp/writer.hpp"

using namespace std;
using namespace libsnark;

// A * B = C;

const int N = 20;

int main()
{
    Writer write("QMP_only_matrix_mul_924.txt");
    Timer timer;
    default_r1cs_gg_ppzksnark_pp::init_public_params();
    typedef libff::Fr<default_r1cs_gg_ppzksnark_pp> FieldT;

    log_debug("QMP");
    libff::enter_block("QMP phrase Start");
    vector<FieldT> A0(N * N, 1);
    vector<FieldT> A1(N * N, 1);
    vector<FieldT> A2(N * N, 1);

    random_generator<FieldT>(A1, N * N);
    sleep(1);
    random_generator<FieldT>(A2, N * N);

    vector<FieldT> A3(N * N, 0);
    multi2(A1, A2, A3, N, N, N);

    log_debug("打印矩阵 A1");
    print_matrix<FieldT>(A1, N, N);
    log_debug("打印矩阵 A2");
    print_matrix<FieldT>(A2, N, N);
    log_debug("打印矩阵 A3");
    print_matrix<FieldT>(A3, N, N);

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
    for (int i = 0; i < cs.A1.size(); ++i)
    {
        g_A_vec[i] = cs.A1[i] * crs.g;
    }

    for (int i = 0; i < cs.A2.size(); ++i)
    {
        g_B_vec[i] = cs.A2[i] * crs.h;
    }
    log_debug("auxiliary A2 generate");

    for (int i = 0; i < cs.A3.size(); ++i)
    {
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
    log_clear();
}
