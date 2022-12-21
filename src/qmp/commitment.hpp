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
commitment_parameter<ppT> commitment_setup(size_t d, size_t l);

template <typename ppT>
commitment<ppT> commitment_generate(polynomial<ppT> &poly, commitment_parameter<ppT> &commitment_parameter, size_t d, size_t l);

template <typename ppT>
bool commitment_verify(commitment<ppT> &C, commitment_parameter<ppT> &commitment_parameter);

template <typename ppT>
void print_Hash_to_file(u<ppT> &u, commitment<ppT> &D, libff::GT<ppT> &U, std::string pathToFile);

std::string print_hash(std::vector<unsigned char> &s1);

std::string hex_to_dec(std::string &s);
template <typename ppT>
prove<ppT> commitment_prove(commitment_parameter<ppT> &commitment_parameter,
                            u<ppT> u,
                            w<ppT> w,
                            int d,
                            int l);

template <typename ppT>
bool commitment_verify_prove(commitment_parameter<ppT> &commitment_parameter,
                             u<ppT> &u,
                             prove<ppT> &pi);