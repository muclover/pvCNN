#ifndef LAYER_H
#define LAYER_H
#include <string>
#include "sealencryptor.h"
#include <ostream>
#include <fstream>

using namespace std;

class Layer
{
public: 
	string name;
	Layer();
    Layer(string layer_name);
	virtual ~Layer();

	string getName(){return name;}
	virtual void printLayerStructure()=0;
	virtual ciphertext3D forward (ciphertext3D input)=0;
	/*If a layer has Plaintext parameters it saves/ loades it in file_name*/
	virtual void savePlaintextParameters(ostream * outfile)=0;
	virtual void loadPlaintextParameters(istream * infile)=0;
	
	/* Compute the last x and y coordinates of a xd*yd image in which is possible to
	apply a filter xf*yf, given a stride [xs,ys].
	Simirarly, in the case of xp*yp pooling the function achieves the same results,
	by substituting them into xf and yf. */
	void computeBoundaries(int xd, int yd, int xs, int ys, int xf, int yf, int* xl, int* yl);

};

#endif