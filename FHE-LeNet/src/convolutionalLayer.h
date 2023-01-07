#ifndef CONVOLUTIONAL_LAYER_H
#define CONVOLUTIONAL_LAYER_H
#include <vector>
#include <string>
#include"seal/seal.h"
#include "layer.h"
#include "sealencryptor.h"

using namespace seal;
class ConvolutionalLayer : public Layer
{
	/* This class implements a three-dimensional convolution on a Ciphertext
	three-dimensional image, characterized by xd, yd and zd.
	The convolution is applied with a stride defined by xs and ys (the z dimension
	has a fixed stride of 1) and with a xf*yf*zf*nf filter, passed as a further
	parameter, with the number f of the filter.
	In addition the class requires the bias for the f-th filter and the average
	image, that is an array with one element per channel.
	The output is the resulting convolved ciphertext image, along with the output dimensions
	xo and yo, that has to be computed with the above function.
	d=dim
	s=stride
	f=filter;
	o=output;
	Note: it is applied only one filter, thus the z-dimension is 1. */

public:
	int xd,yd,zd,xs,ys,xf,yf,nf,th_count;
	int xo,yo,zo;
	plaintext4D filters; //nf,zd,xf,yf
	vector<Plaintext> biases;
	bool filters_already_ntt;
	ConvolutionalLayer(string name,int xd,int yd,int zd,int xs,int ys,int xf,int yf,int nf, int th_count,plaintext4D & filters, vector<Plaintext> & biases);
	ConvolutionalLayer(string name,int xd,int yd,int zd,int xs,int ys,int xf,int yf,int nf,int th_count,istream * infile);
	~ConvolutionalLayer();

	ciphertext3D forward (ciphertext3D input);
	plaintext3D getKernel(int kernel_index);
	Plaintext getBias(int bias_index);
	ciphertext2D convolution3d(ciphertext3D image, plaintext3D kernel,Plaintext bias);
	void transform_input_to_ntt(ciphertext3D &input);
	void transform_kernel_to_ntt(int kernel_index);
	void savePlaintextParameters(ostream * outfile);
	void loadPlaintextParameters(istream * infile);
	void printLayerStructure();
};


#endif