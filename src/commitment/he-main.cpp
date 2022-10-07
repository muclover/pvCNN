#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <bitset>

#include "commitment.hpp"

#include <libff/algebra/fields/field_utils.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark_params.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include <libsnark/gadgetlib1/protoboard.hpp>

#include <seal/seal.h>

// #include "picosha2.h"   //SHA256 lib

using namespace std;
using namespace libsnark;
using namespace seal;
using namespace libff;

template <typename ppT>
void test_hash()
{
    string s = " The first test";
    string hash = picosha2::hash256_hex_string(s);
    cout << hash << endl;

    string f = " The first test";

    picosha2::hash256_one_by_one hasher;
    hasher.init();

    hasher.process(f.begin(), f.end());

    hasher.finish();
    picosha2::get_hash_hex_string(hasher, hash);
    cout << hash << endl;

    // libff::G1<default_r1cs_ppzksnark_pp> g1 = libff::G1<default_r1cs_ppzksnark_pp>::one();
    // libff::G2<default_r1cs_ppzksnark_pp> g2 = libff::G2<default_r1cs_ppzksnark_pp>::one();

    // libff::G1_precomp<default_r1cs_ppzksnark_pp> g1_precomp=default_r1cs_ppzksnark_pp::precompute_G1(g1);
    // libff::G2_precomp<default_r1cs_ppzksnark_pp> g2_precomp=default_r1cs_ppzksnark_pp::precompute_G2(g2);

    // libff::Fqk<default_r1cs_ppzksnark_pp> kc_A_1 = default_r1cs_ppzksnark_pp::miller_loop(g1_precomp, g2_precomp);

    // print_Hash_to_file<default_r1cs_ppzksnark_pp>(u,commit_Q,kc_A_1,"test.txt");
    // print_Hash_to_file<default_r1cs_ppzksnark_pp>(u,commit_Q,kc_A_1,"test1.txt");

    // ifstream f1("test.txt", std::ios::binary);
    // ifstream f2("test.txt", std::ios::binary);
    // vector<unsigned char> s1(picosha2::k_digest_size);
    // vector<unsigned char> s2(picosha2::k_digest_size);

    // picosha2::hash256(f1, s1.begin(), s1.end());

    // for (int i=0; i< s1.size();i++)
    // {
    //     printf("%.2x ",s1[i]);
    // }
    // cout <<endl;
    // cout << s1.size()<<endl;

    // picosha2::hash256(f2, s2.begin(), s2.end());

    // for (int i=0; i< s2.size();i++)
    // {
    //     printf("%.2x ",s2[i]);
    // }
    // cout <<endl;
    // cout << s2.size()<<endl;
}
void saveCiphertext(Ciphertext encrypted, string filename)
{
    ofstream ct;
    ct.open(filename, ios::binary);
    encrypted.save(ct);
};

// int main()
// {
//     typedef libff::Fr<default_r1cs_ppzksnark_pp> FieldT;    //typedef libff::default_ec_pp default_r1cs_ppzksnark_pp; in default_type directory files.
//     default_r1cs_ppzksnark_pp::init_public_params();

//     //d---X l---Y
//     size_t d=200,l=300;

//     FieldT k = 2;
//     // int k =2;

//     polynomial<default_r1cs_ppzksnark_pp> P;
//     for (size_t i=0;i<=d;i++)
//     {
//         P.coeffs.emplace_back(i);
//         for(size_t j=1;j<=l;j++)
//         {
//             P.coeffs.emplace_back(0);
//         }
//     }

//     cout <<"********P**********"<<endl;
//     print_polynomial(P);

//     //Q(Y)=P(k,Y)
//     polynomial<default_r1cs_ppzksnark_pp> Q;
//     FieldT num=FieldT::zero();
//     for(size_t j=0;j<=l;j++)
//     {
//         for(size_t i=0;i<=d;i++)
//         {
//             num += (P.coeffs[i*(l+1)+j]*(k^i));
//         }
//         Q.coeffs.emplace_back(num);
//         num=FieldT::zero();
//     }
//     for(size_t i=l+1;i<P.coeffs.size();i++)
//     {
//         Q.coeffs.emplace_back(0);
//     }
//     cout <<"********Q**********"<<endl;
//     print_polynomial(Q);

//     num=FieldT::zero();
//     //W=(P-Q)/(X-k)
//     polynomial<default_r1cs_ppzksnark_pp> W;
//     for(size_t t=0;t<=d-1;t++)
//     {
//         for(size_t j=0;j<=l;j++)
//         {
//             for(size_t i=1;i<=d-t;i++)
//             {
//                num += P.coeffs[(i+t)*(l+1)+j] * (k^(i-1));
//             }
//             W.coeffs.emplace_back(num);
//             num =FieldT::zero();
//         }

//     }

//     for(size_t i=0;i<=l;i++)
//     {
//         W.coeffs.emplace_back(0);
//     }
//     cout <<"********W**********"<<endl;
//     print_polynomial(W);

//     commitment_parameter<default_r1cs_ppzksnark_pp> parameter;
//     parameter = commitment_setup<default_r1cs_ppzksnark_pp>(d,l);

//     cout <<"********Generate commitment for P**********"<<endl;

//     commitment<default_r1cs_ppzksnark_pp> commit_P;
//     commit_P = commitment_generate<default_r1cs_ppzksnark_pp>(P, parameter,d,l);

//     bool result = commitment_verify<default_r1cs_ppzksnark_pp>(commit_P , parameter);
//     cout << result <<endl;

//     cout <<"*******************************************"<<endl;

//     cout <<"********Generate commitment for Q**********"<<endl;
//     commitment<default_r1cs_ppzksnark_pp> commit_Q;
//     commit_Q = commitment_generate<default_r1cs_ppzksnark_pp>(Q, parameter,d,l);

//     result = commitment_verify<default_r1cs_ppzksnark_pp>(commit_Q , parameter);
//     cout << result <<endl;
//     cout <<"*******************************************"<<endl;
//     FieldT a=FieldT::one()*2;//Field a=Field::one()*k; when k is int type

//     u<default_r1cs_ppzksnark_pp> u(commit_P.c,commit_P.c_hat,commit_Q.c,commit_Q.c_hat,a);
//     w<default_r1cs_ppzksnark_pp> w(P,Q,commit_P.rho,commit_Q.rho);
//     cout <<"*****************===***********************"<<endl;
//     cout <<"*****************===***********************"<<endl;
//     cout <<"*****************===***********************"<<endl;
//     cout <<"*****************===***********************"<<endl;
//     cout <<"*****************===***********************"<<endl;

//     cout << w.P.coeffs.size()<<endl;
//     cout << w.Q.coeffs.size()<<endl;
//     cout <<"*******************************************"<<endl;

//     prove<default_r1cs_ppzksnark_pp> pi=commitment_prove<default_r1cs_ppzksnark_pp>(parameter,u,w,d,l);

//     result = commitment_verify_prove(parameter, u, pi);
//     cout << result <<endl;

// }

/*
Helper function: Prints the parameters in a SEALContext.
*/
inline void print_parameters(const seal::SEALContext &context)
{
    auto &context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

/*
    Helper Function: print matrix
    Notice that it's for test, only suitable for low dimension matrix
*/
template <typename T>
void printf_matrix(vector<vector<T>> v, size_t col)
{
    cout << "**************************" << endl;
    cout << "Matrix Print:" << endl;
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < col; j++)
        {
            cout << " " << v[i][j];
        }
        cout << endl;
    }
    cout << "Matrix Print end!" << endl;
    cout << "**************************" << endl;
}
template <typename T>
void printf_vector(vector<T> v)
{
    cout << "vector :";
    for (int i = 0; i < v.size(); i++)
    {
        cout << " " << v[i];
    }
    cout << endl;
    cout << "vector size: " << v.size();
    cout << endl;
}
void example_mvproduct()
{
    // TODO implement

    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    EncryptionParameters parms(scheme_type::ckks);

    vector<int> moduli(4, 40);
    moduli[0] = 50;
    moduli[moduli.size() - 1] = 59;

    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, moduli));

    double scale = pow(2.0, 40);
    SEALContext context(parms);

    print_parameters(context);

    cout << "Generating keys...";
    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    // Generate the galois keys
    GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);

    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << endl;
    cout << "...done " << endl;

    // generate dim*dim matrix and weight * weight matrix
    // and flat those matrix into vector
    int dim = 50;
    int weight = 5;
    int block = dim - weight + 1;
    int block_size = dim * weight;
    int matrix_size = block * block_size;

    vector<vector<double>> M_dim(dim); // matrix
    vector<double> V_dim;              // size : matrix_size
    int d = 1;
    for (int i = 0; i < M_dim.size(); i++)
    {
        M_dim[i].resize(dim);
        for (int j = 0; j < dim; j++)
        {
            M_dim[i][j] = d++;
        }
    }
    // printf_matrix<double>(M_dim,dim);

    for (int i = 0; i < dim - weight + 1; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            for (int k = 0; k < weight; k++)
            {
                V_dim.push_back(M_dim[i + k][j]);
            }
        }
    }
    // printf_vector<double>(V_dim);

    d = 1;
    vector<vector<double>> M_weight(weight); // weight matrix
    vector<double> V_weight;                 // size : weight*weight
    for (int i = 0; i < M_weight.size(); i++)
    {
        M_weight[i].resize(weight);
        for (int j = 0; j < weight; j++)
        {
            M_weight[i][j] = d++;
        }
    }

    for (int j = 0; j < weight; j++)
    {
        for (int i = 0; i < weight; i++)
        {
            V_weight.push_back(M_weight[i][j]);
        }
    }
    printf_matrix<double>(M_weight, weight);

    V_weight.insert(V_weight.end(), block_size - V_weight.size(), 0);
    // printf_vector<double>(V_weight);

    // generate  Convolution matrix
    vector<vector<double>> M_conv;
    vector<double> tmp(matrix_size - V_weight.size(), 0);

    for (int i = 0; i < block; i++)
    {
        tmp.insert(tmp.begin() + block_size * i, V_weight.begin(), V_weight.end());
        for (int j = 0; j < block; j++)
        {
            M_conv.push_back(tmp);
            rotate(tmp.begin(), tmp.end() - weight, tmp.end());
        }
        tmp.clear();
        tmp.insert(tmp.end(), matrix_size - V_weight.size(), 0);
    }
    for (int i = 0; i < matrix_size - block * block; i++)
    {
        M_conv.push_back(vector<double>(matrix_size, 1));
    }
    // printf_matrix<double>(M_conv,matrix_size);

    // plaintext computation
    vector<double> Mv(matrix_size, 0);
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            Mv[i] += M_conv[i][j] * V_dim[j];
        }
    }
    // printf_vector(Mv);
    cout << matrix_size << endl;
    // Encode the diagonals
    vector<Plaintext> ptxt_diag(matrix_size);
    Plaintext ptxt_vec;
    for (int i = 0; i < matrix_size; i++)
    {
        vector<double> diag(matrix_size);
        for (int j = 0; j < matrix_size; j++)
        {
            diag[j] = M_conv[j][(j + i) % matrix_size];
        }
        encoder.encode(diag, scale, ptxt_diag[i]);
    }

    // vector<double> vrep(encoder.slot_count()-V_dim.size());
    // vrep.insert(vrep.begin(),V_dim.begin(),V_dim.end());
    // printf_vector(vrep);
    // encoder.encode(vrep, scale, ptxt_vec);

    // Ciphertext ctv;
    // encryptor.encrypt(ptxt_vec, ctv);

    // // Now: perform the multiplication
    // Ciphertext temp;
    // Ciphertext enc_result;
    // for (int i =0; i < matrix_size ; i++){
    //         // rotate
    //         evaluator.rotate_vector(ctv, i, galois_keys, temp);
    //         // multiply
    //         evaluator.multiply_plain_inplace(temp, ptxt_diag[i]);
    //         if (i == 0){
    //             enc_result = temp;
    //         }else{
    //             evaluator.add_inplace(enc_result, temp);
    //         }
    // }

    // Plaintext plain_result;
    // vector<double> result;

    // decryptor.decrypt(enc_result, plain_result);
    // encoder.decode(plain_result, result);

    // for (int i = 0; i < block*block; i++){
    //     cout << "actual: " << result[i] << ", expected: " << Mv[i] << endl;
    // }
}

Ciphertext encrypt_my()
{
    // Y= w * x +b
    const int w = 5;
    const int b = 8;

    int x = 4;

    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 1024;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(1024);

    SEALContext context(parms);
    print_parameters(context);

    cout << "Parameter validation (success): " << context.parameter_error_message() << endl;

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    Plaintext w_plain(to_string(w));
    Plaintext b_plain(to_string(b));
    cout << "Express w = " + to_string(w) + " as a plaintext polynomial 0x" + w_plain.to_string() + "." << endl;
    cout << "Express b = " + to_string(b) + " as a plaintext polynomial 0x" + b_plain.to_string() + "." << endl;

    Ciphertext x_encrypted;
    Plaintext x_plain(to_string(x));
    encryptor.encrypt(x_plain, x_encrypted);
    cout << "    + size of freshly encrypted x: " << x_encrypted.size() << endl; // poly number

    cout << "    + noise budget in freshly encrypted x: " << decryptor.invariant_noise_budget(x_encrypted) << " bits"
         << endl;
    Plaintext x_decrypted;
    cout << "    + decryption of x_encrypted: ";
    decryptor.decrypt(x_encrypted, x_decrypted);
    cout << "0x" << x_decrypted.to_string() << " ...... Correct." << endl;
    cout << "Compute w_mul_x_plus_b (w*x+b)." << endl;
    Ciphertext w_mul_x_plus_b = x_encrypted;
    evaluator.multiply_plain_inplace(w_mul_x_plus_b, w_plain);
    evaluator.add_plain_inplace(w_mul_x_plus_b, b_plain);

    Plaintext decrypted_result;
    decryptor.decrypt(w_mul_x_plus_b, decrypted_result);
    cout << "0x" << decrypted_result.to_string() << " ...... Correct." << endl;
    cout << "    + noise budget in freshly encrypted x: " << decryptor.invariant_noise_budget(w_mul_x_plus_b) << " bits"
         << endl;

    auto data = w_mul_x_plus_b.data(0);

    cout << x_encrypted.poly_modulus_degree() << endl;
    cout << x_encrypted.coeff_modulus_size() << endl;
    cout << x_encrypted.size() << endl;

    for (size_t i = 0; i < x_encrypted.poly_modulus_degree(); i++)
    {
        cout << i << "++" << *(data) << endl;
        data++;
    }
    saveCiphertext(x_encrypted, "cipher.txt");
    return w_mul_x_plus_b;
}

template <typename ppT>
void commit_test()
{

    typedef libff::Fr<ppT> FieldT;

    Ciphertext mu = encrypt_my();

    size_t d = mu.poly_modulus_degree();
    size_t l = 10;
    FieldT k = 2;
    FieldT temp = FieldT::zero();
    auto data = mu.data(0);

    polynomial<ppT> P;
    libff::enter_block("Compute the polynomial P");
    for (size_t i = 0; i <= d; i++)
    {
        temp = *(data + i);
        P.coeffs.emplace_back(temp);
        for (size_t j = 1; j <= l; j++)
        {
            P.coeffs.emplace_back(0);
        }
    }
    libff::leave_block("Compute answer to P", false);

    // print_polynomial(P);
    cout << "The size of P " << P.coeffs.size() << endl;
    libff::enter_block("Compute the polynomial Q");

    polynomial<ppT> Q;
    FieldT num = FieldT::zero();

    for (size_t j = 0; j <= l; j++)
    {
        for (size_t i = 0; i <= d; i++)
        {
            num += (P.coeffs[i * (l + 1) + j] * (k ^ i));
        }
        Q.coeffs.emplace_back(num);
        num = FieldT::zero();
    }
    for (size_t i = l + 1; i < P.coeffs.size(); i++)
    {
        Q.coeffs.emplace_back(0);
    }
    libff::leave_block("Compute answer to Q", false);

    cout << "********Q**********" << endl;
    // print_polynomial(Q);

    cout << "The size of Q " << Q.coeffs.size() << endl;

    commitment_parameter<ppT> parameter;
    parameter = commitment_setup<ppT>(d, l);

    cout << "********Generate commitment for P**********" << endl;

    //产生承诺
    commitment<ppT> commit_P;
    commit_P = commitment_generate<ppT>(P, parameter, d, l);
    //验证承诺
    bool result = commitment_verify<ppT>(commit_P, parameter);
    cout << result << endl;

    cout << "*******************************************" << endl;

    cout << "********Generate commitment for Q**********" << endl;
    commitment<ppT> commit_Q;
    commit_Q = commitment_generate<ppT>(Q, parameter, d, l);
    result = commitment_verify<ppT>(commit_Q, parameter);
    cout << result << endl;
    cout << "*******************************************" << endl;

    u<ppT> u(commit_P.c, commit_P.c_hat, commit_Q.c, commit_Q.c_hat, k);
    w<ppT> w(P, Q, commit_P.rho, commit_Q.rho);

    prove<ppT> pi = commitment_prove<ppT>(parameter, u, w, d, l);

    result = commitment_verify_prove(parameter, u, pi);
    cout << result << endl;
}

void zk_test()
{
    typedef libff::Fr<default_r1cs_gg_ppzksnark_pp> FieldT; // typedef libff::default_ec_pp default_r1cs_gg_ppzksnark_pp; in default_type directory files.

    // Create protoboard
    protoboard<FieldT> pb;

    // Define variables

    pb_variable<FieldT> a;
    pb_variable<FieldT> b;
    pb_variable<FieldT> c;
    pb_variable<FieldT> y1;
    pb_variable<FieldT> y2;
    pb_variable<FieldT> y;

    // Allocate variables to protoboard
    // The strings (like "x") are only for debugging purposes
    y.allocate(pb, "y");
    a.allocate(pb, "a");
    b.allocate(pb, "b");
    c.allocate(pb, "c");
    y1.allocate(pb, "y1");
    y2.allocate(pb, "y2");

    // This sets up the protoboard variables
    // so that the first one (out) represents the public
    // input and the rest is private input
    pb.set_input_sizes(1);

    // Add R1CS constraints to protoboard

    // a*b = y1
    pb.add_r1cs_constraint(r1cs_constraint<FieldT>(a, b, y1));

    // c * 1 = y2
    pb.add_r1cs_constraint(r1cs_constraint<FieldT>(c, 1, y2));

    // (y1+y2) * (y1+y2) = y
    pb.add_r1cs_constraint(r1cs_constraint<FieldT>(y1 + y2, y1 + y2, y));

    // Add witness values

    pb.val(a) = 3;
    pb.val(b) = 2;
    pb.val(c) = 4;
    pb.val(y) = 100;
    pb.val(y1) = 6;
    pb.val(y2) = 4;

    const r1cs_constraint_system<FieldT> constraint_system = pb.get_constraint_system();

    const r1cs_gg_ppzksnark_keypair<default_r1cs_gg_ppzksnark_pp> keypair = r1cs_gg_ppzksnark_generator<default_r1cs_gg_ppzksnark_pp>(constraint_system);

    const r1cs_gg_ppzksnark_proof<default_r1cs_gg_ppzksnark_pp> proof = r1cs_gg_ppzksnark_prover<default_r1cs_gg_ppzksnark_pp>(keypair.pk, pb.primary_input(), pb.auxiliary_input());

    bool verified = r1cs_gg_ppzksnark_verifier_strong_IC<default_r1cs_gg_ppzksnark_pp>(keypair.vk, pb.primary_input(), proof);

    cout << "Number of R1CS constraints: " << constraint_system.num_constraints() << endl;
    cout << "Primary (public) input: " << pb.primary_input() << endl;
    cout << "Auxiliary (private) input: " << pb.auxiliary_input() << endl;
    cout << "Verification status: " << verified << endl;
}

int main()
{
    default_r1cs_gg_ppzksnark_pp::init_public_params();
    typedef libff::Fr<default_r1cs_gg_ppzksnark_pp> FieldT; // typedef libff::default_ec_pp default_r1cs_gg_ppzksnark_pp; in default_type directory files.

    cout << "test" << endl;

    // encrypt_my();
    // example_mvproduct();

    commit_test<default_r1cs_gg_ppzksnark_pp>();
    // zk_test();

    return 0;
}
