#ifndef NETWORK
#define NETWORK

#include "sealencryptor.h"
#include "layer.h"
#include <vector>
#include <exception>

using namespace std;

class Network
{
public:

	vector<shared_ptr<Layer>> layers;

	Network();
	~Network();

	int getNumLayers(){return layers.size();}
	virtual shared_ptr<Layer> getLayer(int i){ return layers[i];}
	vector<shared_ptr<Layer> > & getLayers(){return layers;}
	void printNetworkStructure();
	ciphertext3D forward (ciphertext3D input);
	
};

class OutOfBudgetException : public exception
{
	public:
		const int last_layer_computed;
		OutOfBudgetException(int last_layer_computed):
		last_layer_computed(last_layer_computed){}
		virtual const char* what() const throw()
  		{

    	return ("OutOfBudgetException at layer "+to_string(last_layer_computed)).c_str();
  		}
};



#endif