#include "cnnBuilder.h"
#include "sealencryptor.h"
#include "convolutionalLayer.h"
#include "fullyConnectedLayer.h"
#include "network.h"
#include <string>
#include <vector>
#include <cmath>
#include <ostream>
#include <fstream>

using namespace std;

	CnnBuilder::CnnBuilder(string plain_model_path){
		
		ldata.setFileName(plain_model_path);

	}

	vector<float> CnnBuilder::getPretrained(string var_name){
		ldata.setVarName(var_name);
   		return ldata.getData();
	}

	ConvolutionalLayer * CnnBuilder::buildConvolutionalLayer(string name,int xd,int yd,int zd,int xs,int ys,int xf,int yf,int nf, int th_count,istream * infile){
		if(infile!=NULL)
			return new ConvolutionalLayer(name,xd,yd,zd,xs,ys,xf,yf,nf,th_count,infile);

		int n,z,i,j,w;
		vector<float> weights,biases;
		plaintext4D encoded_weights(nf, plaintext3D(zd, plaintext2D (xf, vector<Plaintext> (yf ) )));
		vector<Plaintext> encoded_biases(nf);

		weights=getPretrained(name+".weight");
		biases=getPretrained(name+".bias");
		w=0;
		for(n=0;n<nf;n++){
			for(z=0;z<zd;z++)
				for(i=0;i<xf;i++)
					for(j=0;j<yf;j++){
						encoded_weights[n][z][i][j]=fraencoder->encode( weights[w]);
						w++;
					}
			encoded_biases[n]=fraencoder->encode(biases[n]);
		}

		return new ConvolutionalLayer(name,xd,yd,zd,xs,ys,xf,yf,nf,th_count,encoded_weights,encoded_biases);

	}


	FullyConnectedLayer * CnnBuilder::buildFullyConnectedLayer(string name, int in_dim, int out_dim, int th_count,istream * infile){
		if(infile!=NULL)
			return new FullyConnectedLayer(name,in_dim,out_dim,th_count,infile);

		int i,j,w;
		vector<float> weights,biases;
		vector<Plaintext> encoded_biases(out_dim);
		plaintext2D encoded_weights(out_dim,vector<Plaintext> (in_dim));

		weights=getPretrained(name+".weight");
		biases=getPretrained(name+".bias");
		w=0;
		for(i=0;i<out_dim;i++){
			for(j=0;j<in_dim;j++){
				encoded_weights[i][j]=fraencoder->encode( weights[w] );
				//cout<<n<<" "<<z<<" "<<i<<" "<<j<<" -> "<<weights[w]<<endl;
				w++;
			}
			encoded_biases[i]=fraencoder->encode(biases[i]);
		}
	return new FullyConnectedLayer(name,in_dim,out_dim,th_count,encoded_weights,encoded_biases);

	}


	PoolingLayer * CnnBuilder::buildPoolingLayer(string name,int xd,int yd, int zd, int xs, int ys,int xf, int yf){
		return new PoolingLayer(name,xd,yd,zd,xs,ys,xf,yf);
	}
	AvgPoolingLayer * CnnBuilder::buildAvgPoolingLayer(string name,int xd,int yd, int zd, int xs, int ys,int xf, int yf){
		return new AvgPoolingLayer(name,xd,yd,zd,xs,ys,xf,yf);
	}
	SquareLayer *CnnBuilder::buildSquareLayer(string name, int th_count){
		return new SquareLayer(name,th_count);
	}

	Network CnnBuilder::buildNetwork(string file_name){
		int th_count=40,th_count2=50,th_tiny=32,th_tiny2=42;
		Network net;
		ifstream *infile=NULL;
		if(file_name!=""){
			infile = new ifstream(file_name, ifstream::binary);
		}
		
		//----Lenet.h5---------
        ConvolutionalLayer *conv1 = buildConvolutionalLayer("features.conv1", 32,32,1,1,1,5,5,6,6,infile);
		net.getLayers().push_back(shared_ptr<Layer> (conv1));
		SquareLayer *act1= buildSquareLayer("act1",6);
		net.getLayers().push_back(shared_ptr<Layer> (act1));
        AvgPoolingLayer *pool1 = buildAvgPoolingLayer("pool1",28,28,6,2,2,2,2);
		net.getLayers().push_back(shared_ptr<Layer> (pool1));
        ConvolutionalLayer *conv2= buildConvolutionalLayer("features.conv2",14,14,6,1,1,5,5,16,16,infile);
		net.getLayers().push_back(shared_ptr<Layer> (conv2));
		SquareLayer *act2= buildSquareLayer("act2",16);
		net.getLayers().push_back(shared_ptr<Layer> (act2));
		AvgPoolingLayer *pool2= buildAvgPoolingLayer("pool2",10,10,16,2,2,2,2);
		net.getLayers().push_back(shared_ptr<Layer> (pool2));
		FullyConnectedLayer *fc3= buildFullyConnectedLayer("classifier.fc3",5*5*16,120,40,infile);
		net.getLayers().push_back(shared_ptr<Layer> (fc3));
		SquareLayer *act3= buildSquareLayer("act3",40);
		net.getLayers().push_back(shared_ptr<Layer> (act3));		
		FullyConnectedLayer *fc4= buildFullyConnectedLayer("classifier.fc4",120,84,40,infile);
		net.getLayers().push_back(shared_ptr<Layer> (fc4));
		SquareLayer *act4= buildSquareLayer("act4",40);
		net.getLayers().push_back(shared_ptr<Layer> (act4));
		FullyConnectedLayer *fc5= buildFullyConnectedLayer("classifier.fc5",84,10,40,infile);
		net.getLayers().push_back(shared_ptr<Layer> (fc5));	

		if(infile!=NULL){
			infile->close();
			delete infile;
		}

		return net;
	}
	//Precondition: setAndSaveParamters() or initFromKeys() must have been called before
	Network CnnBuilder::buildAndSaveNetwork(string file_name){
		ofstream * outfile = new ofstream(file_name, ofstream::binary);

		Network net=buildNetwork();
		for(int i=0; i<net.getNumLayers();i++){
			//cerr<<i<<endl<<flush;
			net.getLayer(i)->savePlaintextParameters(outfile);
		}

		outfile->close();

		delete outfile;

		return net;

	}

    //Precondition: setAndSaveParamters() or initFromKeys() must have been called before
	//Save network eac layer
	Network CnnBuilder::buildAndSaveLayer(string file_name, int layerIndex){
		ofstream * outfile = new ofstream(file_name, ofstream::binary);
        Network net=buildNetwork();
		net.getLayer(layerIndex)->savePlaintextParameters(outfile);
		outfile->close();
		delete outfile;
		return net;
	}

	CnnBuilder::~CnnBuilder(){};