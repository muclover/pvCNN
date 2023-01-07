#include "network.h"
#include "sealencryptor.h"
#include "seal/seal.h"
#include "layer.h"
#include <cassert>
#include <chrono>
using namespace std;
using namespace seal;


	Network::Network(){}
	
	void Network::printNetworkStructure(){
		for(int i=0; i<layers.size();i++){
			cerr<<"("<<i<<") : ";
			layers[i]->printLayerStructure();
			cout<<endl;
		}
		

	}
	ciphertext3D Network::forward (ciphertext3D input){
		cout<<"Forward!"<<endl;

		int layer_before_reenc=6;		
		// chrono::high_resolution_clock::time_point time_start, time_end;
		// chrono::microseconds time_forward(0);
		// chrono::microseconds time_reenc(0);
        // time_start= chrono::high_resolution_clock::now();
		cout<<"Layers forward!"<<endl;
        //layers.size()
		for(int i=0; i<8;i++){
			if(i==layer_before_reenc){
				//time_start= chrono::high_resolution_clock::now();
				floatCube image=decryptImage(input);
				
				input=encryptImage(image);
				//time_end= chrono::high_resolution_clock::now();
				//time_reenc =chrono::duration_cast<chrono::microseconds>(time_end - time_start);
				//cout<<time_reenc.count()<<",";

			}
			      
			input=layers[i]->forward(input);
			
			}
        // time_end= chrono::high_resolution_clock::now();
		// time_forward=chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        // cout<<"Forward Time: "<<time_forward.count()<<"!"<<endl;
		return input;
	}

	
//Forward for the optimalParametersChooser---------

	// ciphertext3D Network::forward (ciphertext3D input){
	// 	/*Necessary when there is the need to find optimal plain_modulus*/
	// 	/*If noise budget after the computation of the forward on a layer, is zero, than we reencrypt the output of the previous layer
	// 	and repeat the forward for that layer*/
	// 	int max_num_of_reencryptions=1;
	// 	ciphertext3D output;
	// 	//chrono::high_resolution_clock::time_point time_start, time_end;


	// 	for(int i=0; i<layers.size();i++){
	// 		cerr<<i<<" "<<decryptor->invariant_noise_budget(input[0][0][0])<<endl;
	// 		//chrono::microseconds time_forward(0);
	// 		//chrono::microseconds time_reenc(0);
	// 		//time_start= chrono::high_resolution_clock::now();
	// 		output=layers[i]->forward(input);
	// 		//time_end= chrono::high_resolution_clock::now();
	// 		//time_forward+=chrono::duration_cast<chrono::microseconds>(time_end - time_start);
			
	// 		if(decryptor->invariant_noise_budget(output[0][0][0])<=5){
	// 			if(max_num_of_reencryptions>0){
	// 				//time_start= chrono::high_resolution_clock::now();
	// 				floatCube image=decryptImage(input);
	// 				input=encryptImage(image);
	// 				//time_end= chrono::high_resolution_clock::now();
	// 				//time_reenc =chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	// 				//cout<<time_forward.count()<<","<<time_reenc.count()<<",";
	// 				i--;
	// 				max_num_of_reencryptions--;
	// 				cout<<"Reencrypt due to out of noise budget. Still remaining "<<max_num_of_reencryptions<<endl;
	// 				continue;
	// 			}
	// 			/* Used for optimalParametrsChooser: If max num of reencryptions has been reached throw exception telling which is the maximum layer that can be computed correclty*/
	// 			else{
	// 			 	throw OutOfBudgetException(i-1);
	// 			}
	// 		}
	// 		//time_start= chrono::high_resolution_clock::now();
	// 		input=deepCopyImage(output);
	// 		//time_end= chrono::high_resolution_clock::now();
	// 		//time_forward+=chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	// 		//cout<<time_forward.count()<<","<<"0"<<",";

	// 	}
	// 	return output;
	// }

//------------------------------------

	Network::~Network(){}

