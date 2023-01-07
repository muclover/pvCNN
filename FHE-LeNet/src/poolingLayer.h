#ifndef POOLING_LAYER
#define POOLING_LAYER

#include "seal/seal.h"
#include "sealencryptor.h"
#include "layer.h"


using namespace seal;

class PoolingLayer: public Layer
{ /* This class implements a 2d "sum" Pooling on a xd by yd image with a filter("of ones") of size xf by yf, given xs and ys strides. */
public:
	int xd,yd,xs,ys,xf,yf,xo,yo,zo;
	PoolingLayer(string name,int xd,int yd, int zd, int xs, int ys,int xf, int yf);
	~PoolingLayer();

	ciphertext3D forward (ciphertext3D input);
	void printLayerStructure();
	// //No parameters, so no need to implement these methods
	void savePlaintextParameters(ostream * outfile){};
	void loadPlaintextParameters(istream * infile){};
	// /*Simulates pooling computation during forward. The input is not a real chipertext image but it is a representation given by
	// an upper bound on the number of non-zero coefficients in the Plaintext polynomial and the max_abs_value of these coefficients 
	// (=1 if freshly encryption with base=3 or base=2)
	// Input=The sim_input simulates an image and has just one pixel for each channel in input.
	// Computation=  each pixel in the sim_input is copied xf*yf times and then sum up them, because we suppose that all the pixels belonging to the same 
	// channel are the same. Mathematically it is like multipling the pixel xf*yf. The copy is necessary to simulate exactly all the computations on a single pixel.
	// Output= one pixel for each channel in output (len=sim_input.size())*/
	// static vector<ChooserPoly> poolingSimulator(vector<ChooserPoly> & sim_input, int xf, int yf);
	// static ChooserPoly poolingSimulator(ChooserPoly sim_input, int xf, int yf);
	
};





#endif