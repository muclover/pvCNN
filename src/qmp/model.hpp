#pragma once
#include <vector>

// max pool
struct poolKernel
{
    int size;
    int stride;
};

struct convKernel
{
    int ch_out, ch_in, size, stribe, padding;
    convKernel(int ch_out_, int ch_in_, int size_, int stribe_, int padding_)
        : ch_out(ch_out_), ch_in(ch_in_), size(size_), stribe(stribe_), padding(padding_) {}
    convKernel(int ch_o, int ch_i, int size_)
        : ch_out(ch_o), ch_in(ch_i), size(size_), stribe(1), padding(0) {}
};

struct fcKernel
{
public:
    int ch_out, ch_in;
    fcKernel(int ch_o, int ch_i) : ch_out(ch_o), ch_in(ch_i) {}
};

class model
{
public:
    explicit model(std::vector<std::vector<convKernel>> c,
                   std::vector<poolKernel> p,
                   std::vector<fcKernel> f,
                   int n_in,
                   int n_out,
                   int n) : conv(std::move(c)), pool(std::move(p)), fc(std::move(f)),
                            c_in(n_in), c_out(n_out), size(n) {}

    std::vector<std::vector<convKernel>> conv;
    std::vector<poolKernel> pool;
    std::vector<fcKernel> fc;
    int c_in, c_out;
    int size;
    void init();
    void getModel();
    void proof();
};