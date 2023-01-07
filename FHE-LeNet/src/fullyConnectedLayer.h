#ifndef FULLY_CONNECTED
#define FULLY_CONNECTED
#include <vector>
#include <string>
#include "seal/seal.h"
#include "layer.h"
#include "sealencryptor.h"
#include <ostream>
#include <fstream>

using namespace seal;

class FullyConnectedLayer: public Layer
{
public:
	string name;
	//th_count = # of threads to use to split the workload
	int in_dim, out_dim, th_count;
	plaintext2D weights;
	vector<Plaintext> biases;
	bool weights_already_ntt;
	FullyConnectedLayer(string name, int in_dim, int out_dim,int th_count, plaintext2D  & weights, vector<Plaintext> & biases );
	//Weights and biases are loaded from file_name
	FullyConnectedLayer(string name, int in_dim, int out_dim,int th_count, istream * infile);
	~FullyConnectedLayer();

	ciphertext3D forward (ciphertext3D input);
	/*Given a ciphertext3D z,x,y, reshapes it with z=1, x= z*x*y , y=1. Row-by-row and transform the ciphertext in ntt
	If no rashping is needed, it only transforms the ciphertext to ntt*/
	ciphertext3D reshapeInput(ciphertext3D input);
	void transform_input_to_ntt(ciphertext3D &input);
	Plaintext getWeight(int x_index,int y_index);
	Plaintext getBias(int x_index);
	void savePlaintextParameters(ostream * outfile);
	void loadPlaintextParameters(istream * infile);
	void printLayerStructure();
	/*Simulates fully connected computation during forward. The input is not a real chipertext image but it is a representation given by
	an upper bound on the number of non-zero coefficients in the Plaintext polynomial and the max_abs_value of these coefficients 
	(=1 if freshly encryption with base=3 or base=2) 
	Input=The sim_input simulates an image and has just one pixel for each channel in input.
	Compitation= as a normal forward. Reshape if the dimension in input i smaller than the required one (eg 1st fully connected) by copying the same value of a channel
	to obtain the dimension of an input image.
	Output= as in normal fully connected (array of dimension equal to out_dim)*/
	//static vector<ChooserPoly> fullyConnectedSimulator(vector<ChooserPoly> & sim_input, vector<float> & weights, vector<float> & biases);
	//static ChooserPoly fullyConnectedSimulator(ChooserPoly sim_input, int in_dim);	
	
};





#endif