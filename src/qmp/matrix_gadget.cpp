#include "matrix_gadget.hpp"

template <typename ppT>
std::pair<libff::G1<ppT>, libff::Fr<ppT>> pedersen_commitment(libff::Fr<ppT> x)
{
    libff::Fr<ppT> r = libff::Fr<ppT>::random_element(); // random number
    libff::G1<ppT> g = libff::G1<ppT>::one();
    libff::G1<ppT> h = libff::G1<ppT>::random_element();
    // com(x;r)=(g^x)*(h^r)
    libff::G1<ppT> com = x * g + r * h;
    return std::make_pair(com, r);
}

template <typename FieldT>
void ReLU_gadget<FieldT>::generate_r1cs_constraints()
{
    // this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in, 1, x),
    //                              FMT(this->annotation_prefix, "x"));
    // this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(less,1,FieldT::one()),FMT(this->annotation_prefix,"less"));
    // this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(0, 1, max),
    //                              FMT(this->annotation_prefix, "max"));

    comparison->generate_r1cs_constraints();

    this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in, FieldT::one() - less, sym1),
                                 FMT(this->annotation_prefix, "sym1"));
    this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(max, less, sym2),
                                 FMT(this->annotation_prefix, "sym2"));
    this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(sym1 + sym2, 1, out),
                                 FMT(this->annotation_prefix, "out"));
}

template <typename FieldT>
void ReLU_gadget<FieldT>::generate_r1cs_witness()
{
    // in.evaluate(this->pb);
    // out.evaluate(this->pb);
    /* compute result */
    // this->pb.val(out) = this->pb.val(in);
    // this->pb.val(x) = this->pb.val(in);
    // this->pb.val(max) = this->pb.val(0);
    comparison->generate_r1cs_witness();
    this->pb.val(sym1) = this->pb.val(in) * (FieldT::one() - this->pb.val(less));
    this->pb.val(sym2) = this->pb.val(max) * this->pb.val(less);

    this->pb.val(out) = this->pb.val(sym1) + this->pb.val(sym2);
    // this->pb.val(out) = this->pb.val(in) * (this->pb.val(less));
}

// //ReLU max(0,x) comparison_gadget
// template<typename FieldT>
// class ReLU_gadget_vector : public gadget<FieldT> {
// private :
//     int n=5;
// public:
//     // pb_variable<FieldT> in;
//     // pb_variable<FieldT> out;
//     // pb_variable<FieldT> max;
//     // pb_variable<FieldT> less;
//     // pb_variable<FieldT> sym1;
//     // pb_variable<FieldT> sym2;
//     vector<ReLU_gadget<FieldT>> relu_vector;

//     ReLU_gadget_vector(protoboard<FieldT> &pb,
//                       vector<ReLU_gadget<FieldT> >&relu_vector,
//                       const std::string &annotation_prefix=""):
//         gadget<FieldT>(pb, annotation_prefix), relu_vector(relu_vector)
//     {
//     };

//     void generate_r1cs_constraints();
//     void generate_r1cs_witness();
// };
// template<typename FieldT>
// void ReLU_gadget_vector<FieldT>::generate_r1cs_constraints()
// {
//     for(int i=0;i<n;i++)
//         relu_vector[i].generate_r1cs_constraints();
// }

// template<typename FieldT>
// void ReLU_gadget_vector<FieldT>::generate_r1cs_witness()
// {
//     for(int i=0;i<n;i++)
//         relu_vector[i].generate_r1cs_witness();
// }

template <typename FieldT>
void ReLU_gadget_vector<FieldT>::generate_r1cs_constraints()
{

    for (int i = 0; i < n; i++)
    {
        comparison[i].generate_r1cs_constraints();
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less, sym1[i]),
                                     FMT(this->annotation_prefix, "sym1"));
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(max, less, sym2[i]),
                                     FMT(this->annotation_prefix, "sym2"));
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(sym1[i] + sym2[i], 1, out[i]),
                                     FMT(this->annotation_prefix, "out"));
    }
    cout << "generate constraints complete" << endl;
}

template <typename FieldT>
void ReLU_gadget_vector<FieldT>::generate_r1cs_witness()
{
    for (int i = 0; i < n; i++)
    {
        comparison[i].generate_r1cs_witness();
        this->pb.val(sym1[i]) = this->pb.val(in[i]) * (FieldT::one() - this->pb.val(less));
        this->pb.val(sym2[i]) = this->pb.val(max) * this->pb.val(less);
        this->pb.val(out[i]) = this->pb.val(sym1[i]) + this->pb.val(sym2[i]);
    }
    cout << " generate witness complete!" << endl;
}
