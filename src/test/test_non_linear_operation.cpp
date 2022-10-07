/*
    test non-linear operation such ReLU, pooling
*/
#include <iostream>

#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <libff/algebra/fields/field_utils.hpp>

#include "../qmp/utils.hpp"
#include "../qmp/matrix_gadget.hpp"
#include "../qmp/timer.hpp"
#include "../qmp/writer.hpp"

using namespace std;
using namespace libsnark;

int main()
{
    Writer write("ReLU-Pooling-time");
    Timer timer;

    default_r1cs_gg_ppzksnark_pp::init_public_params();
    typedef libff::Fr<default_r1cs_gg_ppzksnark_pp> FieldT;

    // read output from previous layer

    // activite

    // cout << "***********ReLu proof**************"<<endl;

    // typedef libff::Fr<default_r1cs_gg_ppzksnark_pp> FieldT;

    // default_r1cs_gg_ppzksnark_pp::init_public_params();

    // protoboard<FieldT> pb;

    // pb_variable<FieldT> max;
    // // vector<pb_variable<FieldT>> max(con_output.size());
    // vector<pb_variable<FieldT>> in(output.size());
    // vector<pb_variable<FieldT>> out(output.size());
    // for(int i=0;i<in.size();i++){
    // 	in[i].allocate(pb,"in");
    // 	out[i].allocate(pb,"out");
    // 	// max[i].allocate(pb,"max");
    // }
    // max.allocate(pb,"max");
    // // pb.set_input_sizes(in.size);
    // // pb.val(out)=5;
    // // for(int i=0;i<max.size();i++){
    // // 	pb.val(max[i])=zero_point_output;
    // // }
    // pb.val(max)=zero_point_output;
    // cout <<"generate gadget-"<<endl;
    // vector<ReLU_gadget<FieldT> > relu;
    // // ReLU_gadget<FieldT> relu;
    // for(int i=0;i<in.size();i++){
    // 	// relu.push_back( ReLU_gadget<FieldT> (pb,in[i],out[i],max[i],"relu"));
    // 	relu.push_back( ReLU_gadget<FieldT> (pb,in[i],out[i],max,"relu"));
    // 	relu[i].generate_r1cs_constraints();
    // }

    // cout << " generate constraints"<<endl;

    // const r1cs_constraint_system<FieldT> constraint_system=pb.get_constraint_system();
    // libff::enter_block("generate CRS");
    // const r1cs_gg_ppzksnark_keypair<default_r1cs_gg_ppzksnark_pp> keypair = r1cs_gg_ppzksnark_generator<default_r1cs_gg_ppzksnark_pp>(constraint_system);
    // // Add witness values
    // // pb.val(in)=18;
    // libff::leave_block("generate CRS");

    // for(int i=0;i<in.size();i++){
    // 	pb.val(in[i])=output[i];
    // }

    // for(int i=0;i<in.size();i++)
    // 	relu[i].generate_r1cs_witness();

    // libff::enter_block("generate proof");
    // const r1cs_gg_ppzksnark_proof<default_r1cs_gg_ppzksnark_pp> proof = r1cs_gg_ppzksnark_prover<default_r1cs_gg_ppzksnark_pp>(keypair.pk, pb.primary_input(), pb.auxiliary_input());
    // libff::leave_block("generate proof");

    // libff::enter_block("verify proof");
    // bool verified = r1cs_gg_ppzksnark_verifier_strong_IC<default_r1cs_gg_ppzksnark_pp>(keypair.vk, pb.primary_input(), proof);
    // libff::leave_block("verify proof");
    // cout << "Number of R1CS constraints: " << constraint_system.num_constraints() << endl;
    // // cout << "Primary (public) input: " << pb.primary_input() << endl;
    // // cout << "Auxiliary (private) input: " << pb.auxiliary_input() << endl;
    // cout << "Verification status: " << verified << endl;
    // cout << in.size()<<endl;

    // pooling

    cout << "hello" << endl;
    // cout << "***************pooling started**********" << endl;
    // vector<float> pooling_vector;
    // for (int i = 0; i < block; i += 2)
    // {
    //     for (int j = 0; j < block; j += 2)
    //     {
    //         float sum = (actv_output[i * block + j] + actv_output[i * block + j + 1] + actv_output[(i + 1) * block + j] + actv_output[(i + 1) * block + j + 1]) / 4.0;
    //         // cout <<sum<<endl;
    //         // cout << i*block+j <<" "<<i*block+j+1<<" "<<(i+1)*block+j<<" "<<(i+1)*block+j+1<<endl;
    //         pooling_vector.push_back(sum);
    //     }
    // }

    // int pooling_size = sqrt(pooling_vector.size());
    // for (int i = 0; i < pooling_size; i++)
    // {
    //     for (int j = 0; j < pooling_size; j++)
    //     {
    //         cout << pooling_vector[i * pooling_size + j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "******************************************" << endl;

    // vector<int> pooling;
    // for (int i = 0; i < pooling_vector.size(); i++)
    // {
    //     pooling.push_back(pooling_vector[i]);
    // }
    // cout << pooling.size() << endl;

    // for (int i = 0; i < pooling_size; i++)
    // {
    //     for (int j = 0; j < pooling_size; j++)
    //     {
    //         cout << pooling[i * pooling_size + j] << " ";
    //     }
    //     cout << endl;
    // }

    // cout << pooling.size() << endl;

    // cout << "--------------------------------------" << endl;

    // protoboard<FieldT> pb;

    // // vector<pb_variable<FieldT>> in(pooling_vector.size());
    // // vector<pb_variable<FieldT>> out(pooling.size());
    // pb_variable_array<FieldT> out;
    // out.allocate(pb, 12 * 12, "out");

    // pb_variable_array<FieldT> in;
    // in.allocate(pb, 24 * 24, "in");

    // for (int i = 0; i < in.size(); i++)
    // {
    //     pb.val(in[i]) = actv_output[i];
    // }
    // for (int i = 0; i < out.size(); i++)
    // {
    //     pb.val(out[i]) = 4 * pooling_vector[i];
    // }

    // // pb.set_input_sizes(out.size());
    // cout << "generate gadget-" << endl;
    // Pooling_gadget<FieldT> pool(pb, in, out, "pool");

    // cout << " generate constraints" << endl;
    // pool.generate_r1cs_constraints();

    // const r1cs_constraint_system<FieldT> constraint_system = pb.get_constraint_system();

    // const r1cs_gg_ppzksnark_keypair<default_r1cs_gg_ppzksnark_pp> keypair = r1cs_gg_ppzksnark_generator<default_r1cs_gg_ppzksnark_pp>(constraint_system);
    // cout << " generate witness" << endl;

    // pool.generate_r1cs_witness();

    // const r1cs_gg_ppzksnark_proof<default_r1cs_gg_ppzksnark_pp> proof = r1cs_gg_ppzksnark_prover<default_r1cs_gg_ppzksnark_pp>(keypair.pk, pb.primary_input(), pb.auxiliary_input());
    // bool verified = r1cs_gg_ppzksnark_verifier_strong_IC<default_r1cs_gg_ppzksnark_pp>(keypair.vk, pb.primary_input(), proof);
    // cout << "Number of R1CS constraints: " << constraint_system.num_constraints() << endl;
    // // cout << "Primary (public) input: " << pb.primary_input() << endl;
    // // cout << "Auxiliary (private) input: " << pb.auxiliary_input() << endl;
    // cout << "Verification status: " << verified << endl;
}