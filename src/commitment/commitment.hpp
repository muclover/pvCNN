#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>

#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/field_utils.hpp>
#include <libff/common/default_types/ec_pp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include "picosha2.h"

template <typename ppT>
class polynomial
{
public:
    libff::Fr_vector<ppT> coeffs;
};

template <typename ppT>
void print_polynomial(polynomial<ppT> &P)
{
    std::cout << "The size of Poly : " << P.coeffs.size() << std::endl;
    for (size_t i = 0; i < P.coeffs.size(); ++i)
    {
        std::cout << "[" << i << "]" << P.coeffs[i] << std::endl;
    }
}

template <typename ppT>
class commitment_parameter
{
public:
    libff::G1<ppT> h;
    libff::G2<ppT> g_2;

    libff::G1<ppT> h_hat;
    libff::G2<ppT> g_2_hat;

    libff::G1_vector<ppT> g_i_j; // i <=d j<=l
    libff::G1_vector<ppT> g_i_j_hat;

    libff::G1<ppT> h_1;
    libff::G2<ppT> g_2_1;
    commitment_parameter(){};
    commitment_parameter(libff::G1<ppT> &h,
                         libff::G2<ppT> &g_2,
                         libff::G1<ppT> &h_hat,
                         libff::G2<ppT> &g_2_hat,
                         libff::G1_vector<ppT> &&g_i_j,
                         libff::G1_vector<ppT> &&g_i_j_hat,
                         libff::G1<ppT> &h_1,
                         libff::G2<ppT> &g_2_1) : h(h),
                                                  g_2(g_2),
                                                  h_hat(h_hat),
                                                  g_2_hat(g_2_hat),
                                                  g_i_j(std::move(g_i_j)),
                                                  g_i_j_hat(std::move(g_i_j_hat)),
                                                  h_1(h_1),
                                                  g_2_1(g_2_1){};
};

template <typename ppT>
class commitment
{
public:
    libff::G1<ppT> c;
    libff::G1<ppT> c_hat;
    libff::Fr<ppT> rho;
    commitment() = default;
    commitment(const libff::G1<ppT> &c,
               libff::G1<ppT> &c_hat,
               libff::Fr<ppT> &rho) : c(c),
                                      c_hat(c_hat),
                                      rho(rho)
    {
    }
};

template <typename ppT>
class u
{
public:
    libff::G1<ppT> c1;
    libff::G1<ppT> c1_hat;
    libff::G1<ppT> c2;
    libff::G1<ppT> c2_hat;
    libff::Fr<ppT> k;
    u() = default;
    u(const libff::G1<ppT> &c1,
      libff::G1<ppT> &c1_hat,
      libff::G1<ppT> &c2,
      libff::G1<ppT> &c2_hat,
      libff::Fr<ppT> &k) : c1(c1),
                           c1_hat(c1_hat),
                           c2(c2),
                           c2_hat(c2_hat),
                           k(k)
    {
    }
};

template <typename ppT>
class w
{
public:
    // libff::G1_vector<ppT> &&g_i_j_hat,
    // g_i_j(std::move(g_i_j)),

    polynomial<ppT> P;
    polynomial<ppT> Q;
    libff::Fr<ppT> rho1;
    libff::Fr<ppT> rho2;
    w() = default;
    w(const polynomial<ppT> &P,
      polynomial<ppT> &Q,
      libff::Fr<ppT> &rho1,
      libff::Fr<ppT> &rho2) : P(P),
                              Q(Q),
                              rho1(rho1),
                              rho2(rho2)
    {
    }
    // w(polynomial<ppT> &P,polynomial<ppT> &Q, libff::Fr<ppT> &rho1, libff::Fr<ppT> &rho2)
    // {
    //     this.P.coeffs=P.coeffs;
    //     this.Q.coeffs=Q.coeffs;
    //     this.rho1=rho1;
    //     this.rho2=rho2;
    // }
};

template <typename ppT>
class prove
{
public:
    libff::G1<ppT> c;
    libff::G1<ppT> c_hat; // D=(c,c_hat)
    libff::Fr<ppT> e;
    libff::Fr<ppT> delta;
    libff::Fr<ppT> tau;
    prove() = default;
    prove(const libff::G1<ppT> &c,
          libff::G1<ppT> &c_hat,
          libff::Fr<ppT> &e,
          libff::Fr<ppT> &delta,
          libff::Fr<ppT> &tau) : c(c),
                                 c_hat(c_hat),
                                 e(e),
                                 delta(delta),
                                 tau(tau)
    {
    }
};

template <typename ppT>
commitment_parameter<ppT> commitment_setup(size_t d, size_t l)
{
    const libff::Fr<ppT> alpha = libff::Fr<ppT>::random_element(),
                         s = libff::Fr<ppT>::random_element(),
                         t = libff::Fr<ppT>::random_element();

    libff::Fr<ppT> h_fr = libff::Fr<ppT>::random_element();

    libff::G1<ppT> g = libff::G1<ppT>::one(); // g
    libff::G1<ppT> g_hat = alpha * g;         // g' = g^ alpha

    libff::G1<ppT> h = h_fr * libff::G1<ppT>::one(); // h = g^ h_fr
    // libff::G1<ppT> h= libff::G1<ppT>::one();   //h = g^ h_fr

    libff::G1<ppT> h_hat = alpha * h; // h^ alpha

    libff::G2<ppT> g_2 = h_fr * libff::G2<ppT>::one(); // g_2
    // libff::G2<ppT> g_2 =libff::G2<ppT>::one();   //g_2

    libff::G2<ppT> g_2_hat = alpha * g_2; // g_2^alpha

    libff::G1_vector<ppT> g_i_j;
    for (size_t i = 0; i <= d; i++)
    {
        for (size_t j = 0; j <= l; j++)
        {
            g_i_j.emplace_back(((s ^ i) * (t ^ j)) * g);
        }
    }
    std::cout << "The size of g_j_j = " << g_i_j.size() << std::endl;
    libff::G1_vector<ppT> g_i_j_hat;
    for (size_t i = 0; i < g_i_j.size(); i++)
    {
        g_i_j_hat.emplace_back(alpha * g_i_j[i]);
    }

    libff::G2<ppT> g_2_1 = s * g_2;
    libff::G1<ppT> h_1 = s * h;
    commitment_parameter<ppT> parameter = commitment_parameter<ppT>(
        h,
        g_2,
        h_hat,
        g_2_hat,
        std::move(g_i_j),
        std::move(g_i_j_hat),
        h_1,
        g_2_1);
    return parameter;
}

template <typename ppT>
commitment<ppT> commitment_generate(polynomial<ppT> &poly, commitment_parameter<ppT> &commitment_parameter, size_t d, size_t l)
{
    libff::Fr<ppT> rho = libff::Fr<ppT>::random_element();
    libff::G1<ppT> c = libff::G1<ppT>::zero(); // h^rho

    for (size_t i = 0; i < poly.coeffs.size(); i++)
    {
        c = c + poly.coeffs[i] * commitment_parameter.g_i_j[i]; //\prod g_i_j^a_i_j
    }
    c = c + rho * commitment_parameter.h;
    libff::G1<ppT> c_hat = libff::G1<ppT>::zero();

    for (size_t i = 0; i < poly.coeffs.size(); i++)
    {
        c_hat = c_hat + poly.coeffs[i] * commitment_parameter.g_i_j_hat[i]; //\prod g_i_j^a_i_j
    }

    c_hat = c_hat + rho * commitment_parameter.h_hat;
    commitment<ppT> commit = commitment<ppT>(c, c_hat, rho);

    return commit;
}

template <typename ppT>
bool commitment_verify(commitment<ppT> &C, commitment_parameter<ppT> &commitment_parameter)
{
    bool result = true;
    libff::G1_precomp<ppT> c_precomp = ppT::precompute_G1(C.c);
    libff::G1_precomp<ppT> c_hat_precomp = ppT::precompute_G1(C.c_hat);
    libff::G2_precomp<ppT> g_2_precomp = ppT::precompute_G2(commitment_parameter.g_2);
    libff::G2_precomp<ppT> g_2_hat_precomp = ppT::precompute_G2(commitment_parameter.g_2_hat);

    libff::Fqk<ppT> kc_A_1 = ppT::miller_loop(c_precomp, g_2_hat_precomp);
    libff::Fqk<ppT> kc_A_2 = ppT::miller_loop(c_hat_precomp, g_2_precomp);
    libff::GT<ppT> kc_A = ppT::final_exponentiation(kc_A_1 * kc_A_2.unitary_inverse());

    if (kc_A != libff::GT<ppT>::one())
    {
        if (!libff::inhibit_profiling_info)
        {
            libff::print_indent();
            printf("pairing verify incorrect.\n");
        }
        result = false;
    }
    printf("Pairing verify correct!\n");
    return result;
}

// // template<typename ppT>
// // bool commitment_open(commitmeprint_polynomialnt<ppT> &C, polynomial<ppT> &poly,commitment_parameter<ppT> &commitment_parameter)
// // {
// //    libff::Fr<ppT> k = libff::Fr<ppT>::random_element();

// // }

template <typename ppT>
void print_Hash_to_file(u<ppT> &u, commitment<ppT> &D, libff::GT<ppT> &U, std::string pathToFile)
{
    std::ofstream vk_data;
    vk_data.open(pathToFile);

    libff::G1<ppT> c1 = u.c1;
    c1.to_affine_coordinates();
    libff::G1<ppT> c1_hat = u.c1_hat;
    c1_hat.to_affine_coordinates();

    libff::G1<ppT> c2 = u.c2;
    c2.to_affine_coordinates();
    libff::G1<ppT> c2_hat = u.c2_hat;
    c2_hat.to_affine_coordinates();

    libff::Fr<ppT> k = u.k;

    libff::G1<ppT> D_c = D.c;
    D_c.to_affine_coordinates();
    libff::G1<ppT> D_c_hat = D.c_hat;
    D_c_hat.to_affine_coordinates();

    vk_data << c1.X << std::endl;
    vk_data << c1.Y << std::endl;

    vk_data << c1_hat.X << std::endl;
    vk_data << c1_hat.Y << std::endl;

    vk_data << c2.X << std::endl;
    vk_data << c2.Y << std::endl;

    vk_data << c2_hat.X << std::endl;
    vk_data << c2_hat.Y << std::endl;

    vk_data << k << std::endl;

    vk_data << D_c.X << std::endl;
    vk_data << D_c.Y << std::endl;
    vk_data << D_c_hat.X << std::endl;
    vk_data << D_c_hat.Y << std::endl;

    vk_data << U << std::endl;

    vk_data.close();
}

std::string print_hash(std::vector<unsigned char> &s1)
{
    char o[3];
    std::stringstream ss;
    for (size_t i = 0; i < s1.size(); i++)
    {
        sprintf(o, "%.2x", s1[i]);
        ss << o;
    }
    std::string s = ss.str();
    // std::cout << "Hash: " << s << std::endl;
    return s;
}

std::string hex_to_dec(std::string &s)
{
    std::string str;
    for (size_t i = 0; i < s.size(); i = i + 2)
    {
        std::string out(s, i, 2);
        int x = stoi(out, nullptr, 16);
        str += std::to_string(x);
    }
    std::cout << "HEX To Dec : " << str << std::endl;
    return str;
}
template <typename ppT>
prove<ppT> commitment_prove(commitment_parameter<ppT> &commitment_parameter,
                            u<ppT> u,
                            w<ppT> w,
                            int d,
                            int l)
{
    std::cout << "********W**********" << std::endl;
    // W=(P-Q)/(X-k)
    libff::Fr<ppT> num = libff::Fr<ppT>::zero();
    polynomial<ppT> W;
    for (size_t t = 0; t <= d - 1; t++)
    {
        for (size_t j = 0; j <= l; j++)
        {
            for (size_t i = 1; i <= d - t; i++)
            {
                num += w.P.coeffs[(i + t) * (l + 1) + j] * (u.k ^ (i - 1));
            }
            W.coeffs.emplace_back(num);
            num = libff::Fr<ppT>::zero();
        }
    }
    for (size_t i = 0; i <= l; i++)
    {
        W.coeffs.emplace_back(0);
    }
    // print_polynomial(W);
    std::cout << "the size of W =" << W.coeffs.size() << std::endl;
    //(D,w) <--BivCom(W)

    commitment<ppT> commit_W = commitment_generate(W, commitment_parameter, d, l); // commit_W=(D,omega)

    bool result = commitment_verify<ppT>(commit_W, commitment_parameter);
    std::cout << "Verify the commitment of W: " << result << std::endl;

    libff::G1<ppT> ghat = commitment_parameter.h_1 - (u.k * commitment_parameter.h);

    libff::Fr<ppT> x = libff::Fr<ppT>::random_element();
    libff::Fr<ppT> y = libff::Fr<ppT>::random_element();

    libff::G1<ppT> g1 = (x * commitment_parameter.h) + (y * ghat);
    libff::G2<ppT> g2 = commitment_parameter.g_2;
    libff::G1_precomp<ppT> g1_precomp = ppT::precompute_G1(g1);
    libff::G2_precomp<ppT> g2_precomp = ppT::precompute_G2(g2);
    libff::Fqk<ppT> U_1 = ppT::miller_loop(g1_precomp, g2_precomp);
    libff::GT<ppT> U = ppT::final_exponentiation(U_1);

    print_Hash_to_file<ppT>(u, commit_W, U, "hash_prove.txt");

    // e=Hash(u,D,U);

    std::ifstream f1("hash_prove.txt", std::ios::binary);
    std::vector<unsigned char> s1(picosha2::k_digest_size);
    picosha2::hash256(f1, s1.begin(), s1.end());
    std::string s = print_hash(s1);
    std::string str = hex_to_dec(s);

    libff::Fr<ppT> e = libff::Fr<ppT>(str.c_str());
    // std::cout <<"e="<<e<<std::endl;

    // delta=x-(rho'-rho)e mod q
    libff::Fr<ppT> delta = x - (w.rho2 - w.rho1) * e;
    // tau = y - we
    libff::Fr<ppT> tau = y - commit_W.rho * e;

    // pi=(D,e,delta,tau)
    prove<ppT> prove(commit_W.c, commit_W.c_hat, e, delta, tau);
    return prove;
}

template <typename ppT>
bool commitment_verify_prove(commitment_parameter<ppT> &commitment_parameter,
                             u<ppT> &u,
                             prove<ppT> &pi)
{
    libff::Fr<ppT> rho = libff::Fr<ppT>::zero();
    commitment<ppT> commitment_P(u.c1, u.c1_hat, rho);
    commitment<ppT> commitment_Q(u.c2, u.c2_hat, rho);

    bool b1 = commitment_verify<ppT>(commitment_P, commitment_parameter);
    bool b2 = commitment_verify<ppT>(commitment_Q, commitment_parameter);

    commitment<ppT> commitment_W(pi.c, pi.c_hat, rho);

    bool b3 = commitment_verify<ppT>(commitment_W, commitment_parameter);
    // A=e(d,g_1/g^k) * e(c/c',g)^{-1}

    libff::G1<ppT> g1 = pi.c;
    libff::G2<ppT> g2 = (commitment_parameter.g_2_1 - (u.k * commitment_parameter.g_2));
    libff::G1_precomp<ppT> g1_precomp = ppT::precompute_G1(g1);
    libff::G2_precomp<ppT> g2_precomp = ppT::precompute_G2(g2);
    libff::Fqk<ppT> A_1 = ppT::miller_loop(g1_precomp, g2_precomp);
    libff::GT<ppT> B_1 = ppT::final_exponentiation(A_1);

    g1 = (u.c1 - u.c2);
    g2 = commitment_parameter.g_2;
    g1_precomp = ppT::precompute_G1(g1);
    g2_precomp = ppT::precompute_G2(g2);
    libff::Fqk<ppT> A_2 = ppT::miller_loop(g1_precomp, g2_precomp);
    libff::GT<ppT> B_2 = ppT::final_exponentiation(A_2.unitary_inverse());

    libff::GT<ppT> B = ppT::final_exponentiation(A_1 * (A_2.unitary_inverse()));

    // U=e(h^delta*g_hat^tau,g_2)A^e
    g1 = (pi.delta * commitment_parameter.h) + (pi.tau * (commitment_parameter.h_1 - (u.k * commitment_parameter.h)));
    g2 = commitment_parameter.g_2;
    g1_precomp = ppT::precompute_G1(g1);
    g2_precomp = ppT::precompute_G2(g2);
    libff::Fqk<ppT> U_1 = ppT::miller_loop(g1_precomp, g2_precomp);
    libff::Fqk<ppT> A = A_1 * A_2.unitary_inverse();
    libff::GT<ppT> U = ppT::final_exponentiation(U_1);
    U = U * (B ^ (pi.e));

    libff::Fr<ppT> random = libff::Fr<ppT>::zero();
    commitment<ppT> commit_W(pi.c, pi.c_hat, random);
    // b4=(e=Hash(u,D,U))
    print_Hash_to_file<ppT>(u, commit_W, U, "hash_verify.txt");

    std::ifstream f1("hash_prove.txt", std::ios::binary);
    std::ifstream f2("hash_verify.txt", std::ios::binary);

    std::vector<unsigned char> s1(picosha2::k_digest_size);
    std::vector<unsigned char> s2(picosha2::k_digest_size);

    picosha2::hash256(f1, s1.begin(), s1.end());
    std::string ss1 = print_hash(s1);

    picosha2::hash256(f2, s2.begin(), s2.end());
    std::string ss2 = print_hash(s2);
    bool b4 = ss1 == ss2;

    f1.close();
    f2.close();

    return b1 & b2 & b3 & b4;
}