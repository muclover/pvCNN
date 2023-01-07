#ifndef  AVG_POOLING_LAYER
#define AVG_POOLING_LAYER

#include "poolingLayer.h"
#include "seal/seal.h"

class AvgPoolingLayer : public PoolingLayer
{/* This class implements a 2d  avgPooling on a xd by yd image with a filter("of ones") of size xf by yf, given xs and ys strides. */
public:
	Plaintext div_factor;
	AvgPoolingLayer(string name,int xd,int yd, int zd, int xs, int ys,int xf, int yf);
	~AvgPoolingLayer();

	ciphertext3D forward (ciphertext3D input);
	
};




#endif