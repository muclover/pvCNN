#include "avgPoolingLayer.h"
#include "seal/seal.h"
#include "layer.h"
#include "sealencryptor.h"
#include "poolingLayer.h"

using namespace std;
using namespace seal;

AvgPoolingLayer::AvgPoolingLayer(string name,int xd,int yd, int zd, int xs, int ys,int xf, int yf):
	PoolingLayer(name,xd,yd,zd,xs,ys,xf,yf){
		div_factor=fraencoder->encode(1./(xf*yf));
}


ciphertext3D AvgPoolingLayer::forward (ciphertext3D input){
	int i,j,kx,ky,z,xlast,ylast,p;
	vector<Ciphertext> pixels(xf*yf);
	Ciphertext tmp;
	ciphertext3D pooled_result(zo,ciphertext2D(xo,vector<Ciphertext>(yo)));
	//Plaintext div;
	// float div_plain=1./(xf*yf);
	// div=fraencoder->encode(div_plain);

	Layer::computeBoundaries(xd,yd,xs,ys,xf,yf, &xlast, &ylast);

	for(z=0; z<zo; z++){
		for(i=0; i<xlast; i+=xs)
			for(j=0; j<ylast; j+=ys){
				p=0;
				for(kx=0; kx<xf; kx++)
					for(ky=0; ky<yf; ky++){
						pixels[p]=Ciphertext(input[z][i+kx][j+ky]);
						p++;
					}

				evaluator->add_many(pixels,pooled_result[z][i/xs][j/ys]);
				evaluator->multiply_plain(pooled_result[z][i/xs][j/ys],div_factor,MemoryPoolHandle::Global());
				

			}
	}

	return pooled_result;
}

AvgPoolingLayer::~AvgPoolingLayer(){}