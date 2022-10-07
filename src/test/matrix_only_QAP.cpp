#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <fstream>
#include <cmath>
#include <vector>

#include "libff/algebra/fields/field_utils.hpp"
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include "libsnark/gadgetlib1/pb_variable.hpp"
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>

#include "../qmp/utils.hpp"
#include "../qmp/matrix_gadget.hpp"
#include "../qmp/writer.hpp"
#include "../qmp/timer.hpp"

using namespace libsnark;
using namespace std;
const int n = 20;

// vector A,B represent two matrix
// C is the product of A and B, also use vector represent
//
// insert A, B and C into protoboard
int main()
{

    Writer write("Time_of_QAP_Matrix_Mul.txt");

    libff::enter_block("QAP phrase");
    typedef libff::Fr<default_r1cs_gg_ppzksnark_pp> FieldT;

    // Initialize the curve parameters

    default_r1cs_gg_ppzksnark_pp::init_public_params();
    protoboard<FieldT> pb;
    vector<pb_variable_array<FieldT>> M_rows(n), N_cols(n); // matrix
    vector<vector<pb_variable<FieldT>>> U;                  // result matrix
    vector<vector<inner_product_gadget<FieldT>>> inner_products;
    // 1. Allocate variables to protoboard
    for (auto i = 0; i < n; i++)
    {
        U.push_back(vector<pb_variable<FieldT>>(n));
    }

    for (auto r = 0; r < n; r++)
    {
        for (auto c = 0; c < n; c++)
        {
            U[r][c].allocate(pb, "U_elt");
        }
    }

    for (auto i = 0; i < n; i++)
    {
        M_rows[i].allocate(pb, n, "M_row");
        N_cols[i].allocate(pb, n, "N_col");
    }
    pb.set_input_sizes(n * n);
    matrix_product_gadget<FieldT> m(pb, n, M_rows, N_cols, U, "matrix_product_gadget");
    m.generate_r1cs_constraints();

    // add witness
    for (auto i = 0; i < n; i++)
    {
        for (auto j = 0; j < n; j++)
        {
            pb.val(M_rows[i][j]) = rand() % (5) + 1;
            pb.val(N_cols[i][j]) = rand() % (5) + 1;
        }
    }
    m.generate_r1cs_witness();
    const r1cs_constraint_system<FieldT> constraint_system = pb.get_constraint_system();

    Timer timer;
    timer.start();

    const r1cs_gg_ppzksnark_keypair<default_r1cs_gg_ppzksnark_pp>
        keypair = r1cs_gg_ppzksnark_generator<default_r1cs_gg_ppzksnark_pp>(constraint_system);
    timer.stop();
    write.keep("*QAP CRS", to_string(timer.elapse_time()));
    timer.clear();

    timer.start();
    const r1cs_gg_ppzksnark_proof<default_r1cs_gg_ppzksnark_pp> proof = r1cs_gg_ppzksnark_prover<default_r1cs_gg_ppzksnark_pp>(keypair.pk, pb.primary_input(), pb.auxiliary_input());
    timer.stop();
    write.keep("*QAP Prover", to_string(timer.elapse_time()));
    timer.clear();

    timer.start();
    bool verified = r1cs_gg_ppzksnark_verifier_strong_IC<default_r1cs_gg_ppzksnark_pp>(keypair.vk, pb.primary_input(), proof);
    timer.stop();
    write.keep("*QAP Verifier", to_string(timer.elapse_time()));
    timer.clear();

    cout << "Number of R1CS constraints: " << constraint_system.num_constraints() << endl;
    // cout << "Primary (public) input: " << (int)pb.primary_input() << endl;      // print result matrix
    // cout << "Auxiliary (private) input: " << pb.auxiliary_input() << endl; // print primary matrix and immediate value
    cout << "Verification status: " << verified << endl;
    write.write(n);
    write.write_string("-------------------- memory --------------");
    write.write_bit("*QAP PK size in bits\t", keypair.pk.size_in_bits());
    write.write_bit("*QAP VK size in bits\t", keypair.vk.size_in_bits());
    write.write_bit("*QAP Proof size in bits\t", proof.size_in_bits());
    libff::leave_block("QAP phrase");
}
