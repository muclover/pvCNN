#ifndef SQUARE_LAYER
#define SQUARE_LAYER

#include "seal/seal.h"
#include "sealencryptor.h"
#include "layer.h"

class SquareLayer: public Layer
{ /*This class  implements an approximate activation function (the square function) by just calling SEAL square function*/
public:
	int th_count;
	SquareLayer(string name,int th_count);
	SquareLayer();
	~SquareLayer();

	ciphertext3D forward (ciphertext3D input);
	void printLayerStructure();
	//No parameters, so no need to implement these methods
	void savePlaintextParameters(ostream * outfile){};
	void loadPlaintextParameters(istream * infile){};
	/*Simulates square computation during forward. The input is not a real chipertext image but it is a representation given by
	an upper bound on the number of non-zero coefficients in the Plaintext polynomial and the max_abs_value of these coefficients 
	(=1 if freshly encryption with base=3 or base=2)
	Input=The sim_input simulates an image and has just one pixel for each channel in input.
	Computation=  as in forward, but only one computation for each channel
	Output= one pixel for each channel in output (len=sim_input.size())*/
	static ChooserPoly squareSimulator(ChooserPoly  sim_input);
	static vector<ChooserPoly> squareSimulator(vector<ChooserPoly> &  sim_input);
};


#endif