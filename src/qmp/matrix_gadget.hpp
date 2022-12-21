#include <iostream>

#include "libsnark/gadgetlib1/gadget.hpp"
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include "libff/algebra/fields/field_utils.hpp"

using namespace libsnark;
using namespace std;

template <typename FieldT>
class matrix_product_gadget : public gadget<FieldT>
{
private:
    vector<vector<inner_product_gadget<FieldT>>> inner_products;
    int n;

public:
    const vector<pb_variable_array<FieldT>> M_rows;
    const vector<pb_variable_array<FieldT>> N_cols; // matrix
    const vector<vector<pb_variable<FieldT>>> U;    // result matrix
    matrix_product_gadget(protoboard<FieldT> &pb,
                          int n_,
                          const vector<pb_variable_array<FieldT>> &M_rows,
                          const vector<pb_variable_array<FieldT>> &N_cols,
                          const vector<vector<pb_variable<FieldT>>> &U,
                          const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), n(n_), M_rows(M_rows), N_cols(N_cols), U(U)
    {
        for (auto i = 0; i < n; i++)
        {
            inner_products.push_back(vector<inner_product_gadget<FieldT>>());
        }

        for (auto r = 0; r < n; r++)
        {
            for (auto c = 0; c < n; c++)
            {
                // U[r][c].allocate(pb, "U_elt");
                inner_products[r].push_back(
                    inner_product_gadget<FieldT>(
                        pb, M_rows[r], N_cols[c], U[r][c], "inner_product"));
            }
        }
    }

    void generate_r1cs_constraints()
    {
        for (auto r = 0; r < n; r++)
        {
            for (auto c = 0; c < n; c++)
            {
                inner_products[r][c].generate_r1cs_constraints();
            }
        }
    }

    void generate_r1cs_witness()
    {
        for (auto r = 0; r < n; r++)
        {
            for (auto c = 0; c < n; c++)
            {
                inner_products[r][c].generate_r1cs_witness();
            }
        }
    }
};

// ReLU max(0,x) comparison_gadget
template <typename FieldT>
class ReLU_gadget : public gadget<FieldT>
{
private:
    // comparison_gadget<FieldT> comparison;
    // pb_variable<FieldT> x;
    // pb_variable<FieldT> max;
    pb_variable<FieldT> less_or_eq;
    std::shared_ptr<comparison_gadget<FieldT>> comparison;

public:
    pb_variable<FieldT> in;
    pb_variable<FieldT> out;
    pb_variable<FieldT> max;
    pb_variable<FieldT> less;
    pb_variable<FieldT> sym1;
    pb_variable<FieldT> sym2;
    ReLU_gadget(){};
    ReLU_gadget(protoboard<FieldT> &pb,
                const pb_variable<FieldT> &in,
                const pb_variable<FieldT> &out,
                const pb_variable<FieldT> &max,
                const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), in(in), out(out), max(max)
    {
        // in.allocate(pb, FMT(this->annotation_prefix, "in"));
        // out.allocate(pb, FMT(this->annotation_prefix, "out"));

        // x.allocate(pb, FMT(this->annotation_prefix, "x"));
        // max.allocate(pb, FMT(this->annotation_prefix, "max"));
        sym1.allocate(pb, "sym_1");
        sym2.allocate(pb, "sym_2");
        less.allocate(pb, FMT(this->annotation_prefix, "less"));
        less_or_eq.allocate(pb, FMT(this->annotation_prefix, "less_or_eq"));

        comparison.reset(new comparison_gadget<FieldT>(pb, 10, in, max, less, less_or_eq,
                                                       FMT(this->annotation_prefix, "comparison")));
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
class ReLU_gadget_vector : public gadget<FieldT>
{
private:
    pb_variable<FieldT> less_or_eq;
    vector<comparison_gadget<FieldT>> comparison;

public:
    int n = 5;
    vector<pb_variable<FieldT>> in{vector<pb_variable<FieldT>>(n)};
    vector<pb_variable<FieldT>> out{vector<pb_variable<FieldT>>(n)};
    pb_variable<FieldT> max;
    pb_variable<FieldT> less;
    vector<pb_variable<FieldT>> sym1{vector<pb_variable<FieldT>>(n)};
    vector<pb_variable<FieldT>> sym2{vector<pb_variable<FieldT>>(n)};

    ReLU_gadget_vector(protoboard<FieldT> &pb,
                       vector<pb_variable<FieldT>> &in,
                       vector<pb_variable<FieldT>> &out,
                       const pb_variable<FieldT> &max,
                       int n,
                       const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), in(in), out(out), max(max), n(n)
    {
        cout << " generate constructor started" << endl;
        for (int i = 0; i < n; i++)
            sym1[i].allocate(pb, "sym_1");
        for (int i = 0; i < n; i++)
            sym2[i].allocate(pb, "sym_2");
        cout << "1" << endl;

        less.allocate(pb, FMT(this->annotation_prefix, "less"));
        less_or_eq.allocate(pb, FMT(this->annotation_prefix, "less_or_eq"));

        for (int i = 0; i < n; i++)
        {
            comparison.push_back(comparison_gadget<FieldT>(pb, 10, in[i], max, less, less_or_eq, "comparison"));
        }
        cout << "constructor complete!" << endl;
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
class Pooling_gadget : public gadget<FieldT>
{
private:
    int n = 24; // for test

public:
    const pb_linear_combination_array<FieldT> in;
    const pb_variable_array<FieldT> out;
    Pooling_gadget(){};
    Pooling_gadget(protoboard<FieldT> &pb,
                   const pb_linear_combination_array<FieldT> &in,
                   const pb_variable_array<FieldT> &out,
                   const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), in(in), out(out)
    {
        cout << " generate constructor started" << endl;
        n = sqrt(in.size());
        cout << "n=" << n << endl;
    };

    void generate_r1cs_constraints()
    {

        int d = 0;
        for (size_t i = 0; i < n; i += 2)
        {
            for (size_t j = 0; j < n; j += 2)
            {

                this->pb.add_r1cs_constraint(
                    r1cs_constraint<FieldT>(in[i * n + j] + in[i * n + j + 1] + in[(i + 1) * n + j] + in[(i + 1) * n + j + 1], 1, out[d++]),
                    FMT(this->annotation_prefix, " S_%zu", i));
            }
        }
    }

    void generate_r1cs_witness()
    {

        for (size_t i = 0; i < in.size(); i++)
        {
            in[i].evaluate(this->pb);
        }
    }
};

template <typename ppT>
std::pair<libff::G1<ppT>, libff::Fr<ppT>> pedersen_commitment(libff::Fr<ppT> x);