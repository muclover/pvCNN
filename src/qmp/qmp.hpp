#pragma once

#include <libff/algebra/curves/public_params.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#include <iostream>

template <typename ppT>
class qmp_zksnark {
public:
    std::vector<libff::Fr<ppT>> A0;
    std::vector<libff::Fr<ppT>> A1;
    std::vector<libff::Fr<ppT>> A2;
    std::vector<libff::Fr<ppT>> A3;
    qmp_zksnark(){};
    qmp_zksnark(std::vector<libff::Fr<ppT>> &&A0,
                std::vector<libff::Fr<ppT>> &&A1,
                std::vector<libff::Fr<ppT>> &&A2,
                std::vector<libff::Fr<ppT>> &&A3) :
        A0(std::move(A0)),
        A1(std::move(A1)),
        A2(std::move(A2)),
        A3(std::move(A3))
    {}
};

/*************** crs ****************/
template <typename ppT>
class qmp_crs {
public:
    libff::G1<ppT> g;
    libff::G2<ppT> h;
    libff::G1<ppT> alpha_g1;
    libff::G1<ppT> beta_g1;
    libff::G2<ppT> beta_g2;
    libff::G1<ppT> delta_g1;
    libff::G2<ppT> delta_g2;
    libff::G2<ppT> gamma_g2;

    libff::G1_vector<ppT> A_query;
    libff::G1_vector<ppT> B_query;
    libff::G1_vector<ppT> C_query;

    qmp_crs()
    {}
    qmp_crs(libff::G1<ppT> &&g,
            libff::G2<ppT> &&h,
            libff::G1<ppT> &&alpha_g1,
            libff::G1<ppT> &&beta_g1,
            libff::G2<ppT> &&beta_g2,
            libff::G1<ppT> &&delta_g1,
            libff::G2<ppT> &&delta_g2,
            libff::G2<ppT> &&gamma_g2,
            libff::G1_vector<ppT> &&A_query,
            libff::G1_vector<ppT> &&B_query,
            libff::G1_vector<ppT> &&C_query) :
        g(std::move(g)),
        h(std::move(h)),
        alpha_g1(std::move(alpha_g1)),
        beta_g1(std::move(beta_g1)),
        beta_g2(std::move(beta_g2)),
        delta_g1(std::move(delta_g1)),
        delta_g2(std::move(delta_g2)),
        gamma_g2(std::move(gamma_g2)),
        A_query(std::move(A_query)),
        B_query(std::move(B_query)),
        C_query(std::move(C_query))
    {}
    size_t G1_size() const
    {
        return 4 + A_query.size() + B_query.size() + C_query.size();
    }
    size_t G2_size() const
    {
        return 4;
    }
    size_t size_in_bits() const
    {
        return (libff::size_in_bits(A_query) + libff::size_in_bits(B_query) + libff::size_in_bits(C_query) + 4 * libff::G1<ppT>::size_in_bits() + 4 * libff::G2<ppT>::size_in_bits());
    }
    void print_size() const
    {
        libff::print_indent();
        printf("* G1 elements in PK: %zu\n", this->G1_size());
        libff::print_indent();
        printf("* G2 elements in PK: %zu\n", this->G2_size());
        libff::print_indent();
        printf("* CRS size in bits: %zu\n", this->size_in_bits());
    }
};

template <typename ppT>
class auxiliary_input_information {
public:
    libff::G1_vector<ppT> g_A_vec;
    libff::G2_vector<ppT> g_B_vec;
    libff::G1_vector<ppT> g_C_vec;
    auxiliary_input_information(){};
    auxiliary_input_information(libff::G1_vector<ppT> &&g_A_vec,
                                libff::G2_vector<ppT> &&g_B_vec,
                                libff::G1_vector<ppT> &&g_C_vec) :
        g_A_vec(std::move(g_A_vec)),
        g_B_vec(std::move(g_B_vec)),
        g_C_vec(std::move(g_C_vec))
    {}

    size_t G1_size() const
    {
        return g_A_vec.size() + g_C_vec.size();
    }

    size_t G2_size() const
    {
        return g_B_vec.size();
    }
    size_t size_in_bits() const
    {
        return (libff::size_in_bits(g_A_vec) + libff::size_in_bits(g_C_vec) + libff::size_in_bits(g_B_vec));
    }

    void print_size() const
    {
        libff::print_indent();
        printf("* G1 elements in Aux: %zu\n", this->G1_size());
        libff::print_indent();
        printf("* G2 elements in Aux: %zu\n", this->G2_size());
        libff::print_indent();
        printf("* Aux size in bits: %zu\n", this->size_in_bits());
    }
};

/*********************************** Proof ***********************************/
/**
 * A proof for the R1CS GG-ppzkSNARK.
 *
 */
template <typename ppT>
class qmp_zksnark_proof {
public:
    // libff::G1_vector<ppT> g_A;
    libff::G1<ppT> g_A;
    // libff::G2_vector<ppT> g_B;
    libff::G2<ppT> g_B;
    libff::G1_vector<ppT> g_C;

    qmp_zksnark_proof(){};

    qmp_zksnark_proof(libff::G1<ppT> &&g_A,
                      libff::G2<ppT> &&g_B,
                      libff::G1_vector<ppT> &&g_C) :
        g_A(std::move(g_A)),
        g_B(std::move(g_B)),
        g_C(std::move(g_C)){};
    size_t G1_size() const
    {
        return 1 + this->g_C.size();
    }

    size_t G2_size() const
    {
        return 1;
    }

    size_t size_in_bits() const
    {
        return G1_size() * libff::G1<ppT>::size_in_bits() + G2_size() * libff::G2<ppT>::size_in_bits();
    }

    void print_size() const
    {
        libff::print_indent();
        printf("* G1 elements in proof: %zu\n", this->G1_size());
        libff::print_indent();
        printf("* G2 elements in proof: %zu\n", this->G2_size());
        libff::print_indent();
        printf("* Proof size in bits: %zu\n", this->size_in_bits());
    }
};

/***************************** Main algorithms *******************************/

/**
 * A prover algorithm for the R1CS GG-ppzkSNARK.
 */
/***************************generate crs****************************/

template <typename ppT>
qmp_crs<ppT> qmp_zksnark_generator(const qmp_zksnark<ppT> &cs, int N)
{
    libff::enter_block("Generate QMP CRS ");

    libff::G1<ppT> g = libff::G1<ppT>::random_element();

    libff::G2<ppT> h = libff::G2<ppT>::random_element();

    /* Generate secret randomness */
    const libff::Fr<ppT> alpha = libff::Fr<ppT>::random_element();
    const libff::Fr<ppT> beta = libff::Fr<ppT>::random_element();
    const libff::Fr<ppT> gamma = libff::Fr<ppT>::random_element();
    const libff::Fr<ppT> delta = libff::Fr<ppT>::random_element();

    const libff::Fr<ppT> gamma_inverse = gamma.inverse();

    libff::G1<ppT> alpha_g1 = alpha * g;
    libff::G1<ppT> beta_g1 = beta * g;
    libff::G2<ppT> beta_g2 = beta * h;
    libff::G1<ppT> delta_g1 = delta * g;
    libff::G2<ppT> delta_g2 = delta * h;

    libff::G2<ppT> gamma_g2 = gamma * h;

    libff::G1_vector<ppT> A_query(cs.A1.size());
    libff::G1_vector<ppT> B_query(cs.A2.size());
    libff::G1_vector<ppT> C_query(cs.A3.size());

    for (auto i = 0; i < N * N; i++) {
        // A_query.emplace_back((beta * gamma_inverse * cs.A1[i])* g);
        // B_query.emplace_back((alpha * gamma_inverse * cs.A2[i])* g);
        // C_query.emplace_back((gamma_inverse * cs.A3[i])* g);
        A_query[i] = ((beta * gamma_inverse * cs.A1[i]) * g);
        B_query[i] = ((alpha * gamma_inverse * cs.A2[i]) * g);
        C_query[i] = ((gamma_inverse * cs.A3[i]) * g);
    }

    libff::leave_block("Generate CRS is done");

    qmp_crs<ppT> crs = qmp_crs<ppT>(std::move(g),
                                    std::move(h),
                                    std::move(alpha_g1),
                                    std::move(beta_g1),
                                    std::move(beta_g2),
                                    std::move(delta_g1),
                                    std::move(delta_g2),
                                    std::move(gamma_g2),
                                    std::move(A_query),
                                    std::move(B_query),
                                    std::move(C_query));
    return crs;
}

/**
 * A prover algorithm for the R1CS GG-ppzkSNARK.
 */
template <typename ppT>
qmp_zksnark_proof<ppT> qmp_zksnark_prover(const qmp_zksnark<ppT> &cs, const qmp_crs<ppT> &crs)
{
    libff::enter_block("Generate Proof Pi ");

    const libff::Fr<ppT> r = libff::Fr<ppT>::random_element();
    const libff::Fr<ppT> u = libff::Fr<ppT>::random_element();

    libff::G1<ppT> g_A = crs.alpha_g1 + (r * crs.delta_g1);
    libff::G2<ppT> g_B = crs.beta_g2 + (u * crs.delta_g2);

    libff::G1<ppT> g_C_unit1 = u * crs.alpha_g1 + r * u * crs.delta_g1;
    libff::G1<ppT> g_C_unit2 = r * crs.beta_g1;

    libff::G1_vector<ppT> g_C(cs.A1.size());

    for (auto i = 0; i < cs.A1.size(); i++) {
        // g_C.emplace_back(g_C1[i]+g_C2[i]);
        g_C[i] = (g_C_unit1 + u * cs.A1[i] * crs.g + g_C_unit2 + r * cs.A2[i] * crs.g);
    }

    libff::leave_block("Generate Proof is done");

    qmp_zksnark_proof<ppT> proof = qmp_zksnark_proof<ppT>(std::move(g_A),
                                                          std::move(g_B),
                                                          std::move(g_C));

    return proof;
}

template <typename ppT>
bool qmp_zksnark_verifier(const qmp_crs<ppT> &crs, const qmp_zksnark_proof<ppT> &proof, auxiliary_input_information<ppT> &aux_input)
{
    bool result = true;
    libff::GT<ppT> alpha_g1_beta_g2 = ppT::reduced_pairing(crs.alpha_g1, crs.beta_g2);

    const libff::G1_precomp<ppT> proof_g_A_precomp = ppT::precompute_G1(proof.g_A);
    const libff::G2_precomp<ppT> proof_g_B_precomp = ppT::precompute_G2(proof.g_B);
    const libff::G2_precomp<ppT> h_precomp = ppT::precompute_G2(crs.h);
    const libff::Fqk<ppT> QMP1_1 = ppT::miller_loop(proof_g_A_precomp, proof_g_B_precomp);

    std::vector<libff::Fqk<ppT>> QMP1_4_vec(aux_input.g_C_vec.size());

    for (auto i = 0; i < aux_input.g_C_vec.size(); i++) {
        libff::G1_precomp<ppT> proof_g_C_vec_precomp = ppT::precompute_G1(aux_input.g_C_vec[i]);
        // QMP1_4_vec.emplace_back( ppT::miller_loop(proof_g_C_vec_precomp,  h_precomp) );
        QMP1_4_vec[i] = (ppT::miller_loop(proof_g_C_vec_precomp, h_precomp));
    }

    for (auto i = 0; i < proof.g_C.size(); i++) {
        // left
        // e(A,B)e(A,h)^[1 2 2 1]*e(g,B)^[2 1 3 2]*e(g,h)^[2 1 3 2]*[1 2 2 1]
        libff::G1_precomp<ppT> proof_g_A_vec_precomp;
        libff::Fqk<ppT> QMP1_3;

        proof_g_A_vec_precomp = ppT::precompute_G1(aux_input.g_A_vec[i]);
        QMP1_3 = ppT::miller_loop(proof_g_A_vec_precomp, proof_g_B_precomp);

        libff::G2_precomp<ppT> proof_g_B_vec_precomp;
        libff::Fqk<ppT> QMP1_2;
        if (aux_input.g_B_vec[i] == libff::G2<ppT>::zero()) {
            QMP1_2 = libff::Fqk<ppT>::one();
        } else {
            proof_g_B_vec_precomp = ppT::precompute_G2(aux_input.g_B_vec[i]);
            QMP1_2 = ppT::miller_loop(proof_g_A_precomp, proof_g_B_vec_precomp);
        }

        libff::Fqk<ppT> QMP1 = QMP1_1 * QMP1_2 * QMP1_3 * QMP1_4_vec[i];

        // right
        libff::G1_precomp<ppT> proof_g_C_precomp = ppT::precompute_G1(proof.g_C[i]);

        libff::G1_precomp<ppT> acc_precomp = ppT::precompute_G1(crs.A_query[i] + crs.B_query[i] + crs.C_query[i]);

        libff::G2_precomp<ppT> gamma_g2_precomp = ppT::precompute_G2(crs.gamma_g2);
        libff::G2_precomp<ppT> delta_g2_precomp = ppT::precompute_G2(crs.delta_g2);

        libff::Fqk<ppT> QMP2_1 = ppT::miller_loop(acc_precomp, gamma_g2_precomp);

        libff::Fqk<ppT> QMP2_2 = ppT::miller_loop(proof_g_C_precomp, delta_g2_precomp);

        libff::Fqk<ppT> QMP2 = QMP2_1 * QMP2_2;
        libff::GT<ppT> QAP = ppT::final_exponentiation(QMP1 * QMP2.unitary_inverse());

        if (QAP != alpha_g1_beta_g2) {
            if (!libff::inhibit_profiling_info) {
                libff::print_indent();
                printf("QMP divisibility check failed.\n");
            }
            result = false;
            break;
        }
    }

    return result;
}
