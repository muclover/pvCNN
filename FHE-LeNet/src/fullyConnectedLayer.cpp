#include "fullyConnectedLayer.h"
#include "layer.h"
#include "seal/seal.h"
#include "sealencryptor.h"
#include <cassert>
#include <thread>
#include <mutex>
#include <ostream>
#include <fstream>

using namespace std;
using namespace seal;


FullyConnectedLayer::FullyConnectedLayer(string name, int in_dim, int out_dim,int th_count, plaintext2D  & weights, vector<Plaintext> & biases):
	Layer(name),
	in_dim(in_dim),out_dim(out_dim),
	th_count(th_count),
	weights(weights), biases(biases),
	weights_already_ntt(false){
	if(th_count>out_dim)
		th_count=out_dim;	
	else if(th_count<=0)
	 		th_count=1;

	}
FullyConnectedLayer::FullyConnectedLayer(string name, int in_dim, int out_dim,int th_count, istream * infile):
	Layer(name),
	in_dim(in_dim),out_dim(out_dim),
	th_count(th_count),
	weights_already_ntt(false){
		loadPlaintextParameters(infile);
		if(th_count>out_dim)
			th_count=out_dim;	
	else if(th_count<=0)
	 		th_count=1;
}
ciphertext3D FullyConnectedLayer::reshapeInput(ciphertext3D input){
	int x_size=input[0].size(), y_size=input[0][0].size(), z_size=input.size();
	//cout<<z_size<<" "<<x_size<<" "<<y_size<<endl;
	int z,x,y;
	if(z_size!=1 && y_size!=1){
		ciphertext3D reshaped_input(1,ciphertext2D(in_dim,vector<Ciphertext>(1)));
		for (int i = 0; i < in_dim; ++i)
		{
			z=i/(x_size*y_size);
			x=i/y_size - (x_size * z);
			y=i%y_size;
			//cout<<z<<" "<<x<<" "<<y<<endl;
			reshaped_input[0][i][0]= Ciphertext(input[z][x][y]);
		}
		return reshaped_input;
	}

	return input;
}

//Transform input in ntt to speedup the multiply_plain, before to start the computation
void FullyConnectedLayer::transform_input_to_ntt(ciphertext3D &input){
    int from=0,to=0;
    int threads=th_count;
    mutex mtx;

    vector<thread> th_vector;


    auto parallelTransform=[&](ciphertext3D &input,mutex &mtx,int from, int to){


        for(int x=from;x<to;x++){
            evaluator->transform_to_ntt(input[0][x][0]);
            }
    };
    
    if(threads>in_dim)
        threads=in_dim;

    int thread_depth=in_dim/threads;

    for (int i = 0; i < threads; i++){
        from=to;
        if(i<threads-1)
            to+=thread_depth;
        else
            to+=thread_depth + (in_dim%threads);

        th_vector.emplace_back(parallelTransform, ref(input),ref(mtx),from,to);
    }

    for (int i = 0; i < threads; i++){
        th_vector[i].join();
    }


}

Plaintext FullyConnectedLayer::getWeight(int x_index,int y_index){
   return weights[x_index][y_index]; 
 }

Plaintext FullyConnectedLayer::getBias(int x_index){
   return biases[x_index]; 
 }
//Transform weight to ntt to speedup the multiply_plain in forward phase
// void FullyConnectedLayer::transform_weights_to_ntt(){
// 	for(int i=0;i<out_dim;i++)
// 		for(int j=0;j<in_dim;j++)
// 			evaluator->transform_to_ntt(weights[i][j],MemoryPoolHandle::Global());

// }

//Forward implemented with threads
ciphertext3D FullyConnectedLayer::forward(ciphertext3D input){

	int from=0,to=0,thread_rows=0;
	ciphertext3D result(1,ciphertext2D(out_dim,vector<Ciphertext>(1)));
	vector<thread> th_vector;

	//Each thread will work on a portion of computation (ax+b) making matrix product of rows from index "form" to intex "to"
	auto parallelForward=[&](ciphertext3D &input,ciphertext3D &result,int from, int to,const MemoryPoolHandle &pool){
		vector<Ciphertext> tmp(in_dim);
		Ciphertext input_copy(pool);

		for(int i=from; i<to; i++){
			for(int j=0;j<in_dim;j++){
				//Temporary copy the input ntt
				input_copy=input[0][j][0];
				//Performed optimized multiply_plain
				if(!weights_already_ntt){
					evaluator->transform_to_ntt(weights[i][j],pool);
				}
				evaluator->multiply_plain_ntt(input_copy,weights[i][j]);
				//Transform the result back to normal Ciphertext
				evaluator->transform_from_ntt(input_copy);
				//Insert the result in tmp vector
				tmp[j]=Ciphertext(input_copy);
			}
		evaluator->add_plain(tmp[0],biases[i]);
		evaluator->add_many(tmp,result[0][i][0]);
		}

	};
	input=reshapeInput(input);
	transform_input_to_ntt(input);

	//transform_weights_to_ntt();

	thread_rows=out_dim/th_count;
	
	
	for (int i = 0; i < th_count; i++){
		from=to;
    	if(i<th_count-1)
    		to+=thread_rows;
    	else
    		to+=thread_rows + (out_dim%th_count);

      	th_vector.emplace_back(parallelForward, ref(input),ref(result),from,to,MemoryPoolHandle::New());

    }
    for (size_t i = 0; i < th_vector.size(); i++)
    {
        th_vector[i].join();
    }
    weights_already_ntt=true;
    return result;

}
//Forward implemented without threads
/*
ciphertext3D FullyConnectedLayer::forward(ciphertext3D input){
	cout<<"Begin forward"<<endl<<flush;
	Plaintext weight;
	ciphertext3D result(1,ciphertext2D(out_dim,vector<Ciphertext>(1)));
	vector<Ciphertext> tmp(in_dim);

	input=reshapeInput(input);

	//cout<<decryptor->invariant_noise_budget(input[0][6][0])<<endl;
	//cout<<input[0].size()<<endl;
	for(int i=0; i<out_dim; i++){
		for(int j=0;j<in_dim;j++){
			weight=getWeight(i,j);
			//cout<<decryptor->invariant_noise_budget(input[0][j][0])<<endl;
			assert(decryptor->invariant_noise_budget(input[0][j][0])>0);
			evaluator->multiply_plain(input[0][j][0],weight,tmp[j]);
			//cout<<decryptor->invariant_noise_budget(tmp[j])<<endl;
			//cout<<"Multiply"<<i<<endl<<flush;
		}
	//adding bias
	evaluator->add_plain(tmp[0],getBias(i));
	evaluator->add_many(tmp,result[0][i][0]);
	//cout<<decryptor->invariant_noise_budget(result[0][i][0])<<endl;
	}
	return result;
}*/

void FullyConnectedLayer::savePlaintextParameters(ostream * outfile){
		int i,j;
		for(i=0;i<out_dim;i++){
			for(j=0;j<in_dim;j++){
				weights[i][j].save(*outfile);
				outfile->flush();
			}
			biases[i].save(*outfile);
			outfile->flush();
		}
}
//Load and transform weights to ntt
void FullyConnectedLayer::loadPlaintextParameters(istream * infile){		
		int i,j;
		vector<Plaintext> encoded_biases(out_dim);
		plaintext2D encoded_weights(out_dim,vector<Plaintext> (in_dim));

		for(i=0;i<out_dim;i++){
			for(j=0;j<in_dim;j++){
				encoded_weights[i][j].load(*infile);
				//weights[i][j].load(infile);
				}
			encoded_biases[i].load(*infile);
			//biases[i].load(infile);
		}
		weights=encoded_weights;
		biases=encoded_biases;
}

void FullyConnectedLayer::printLayerStructure(){
    cerr<<"Fully connected "<<name<<" : ("<<in_dim<<" -> "<<out_dim<<")"<<"run with "<<th_count<<" threads"<<endl;

}

// vector<ChooserPoly> FullyConnectedLayer::fullyConnectedSimulator(vector<ChooserPoly> & sim_input, vector<float> & weights, vector<float> & biases){
// 	cout<<"fc"<<flush;
// 	int approx=1000;
// 	int out_dim=biases.size();
// 	int in_dim=weights.size()/out_dim;
// 	cout<<"in_dim "<<in_dim<<endl;
// 	vector<ChooserPoly> *sim_input_ptr;
// 	vector<ChooserPoly> tmp_sim(in_dim);
// 	vector<ChooserPoly> sim_out(out_dim);

// 	if(sim_input.size()<in_dim){
// 		sim_input_ptr = new vector<ChooserPoly> (in_dim);
// 		int num_copies=in_dim/sim_input.size();
// 		for(int i=0;i<in_dim;i++){
// 			cout<<i/num_copies<<","<<i<<","<<num_copies<<endl;
// 			(*sim_input_ptr)[i]=ChooserPoly(sim_input[int(i/num_copies)]);	
// 		}
// 	}
// 	else sim_input_ptr=&sim_input;

// 	for(int i=0;i<out_dim;i++){
// 		for(int j=0; j<in_dim;j++){
// 	 		cout<<i<<","<<j<<flush<<endl;
// 	 		cout<<"w "<<j+(i*in_dim)<<endl;
// 	 		int weight=weights[j+(i*in_dim)]*approx;
// 			if(weight==0)
// 	 			tmp_sim[j]=chooser_evaluator->multiply_plain((*sim_input_ptr)[j],31,1);
// 	 		else
// 	 			tmp_sim[j]=chooser_evaluator->multiply_plain((*sim_input_ptr)[j],chooser_encoder->encode(weight));
// 	 		//tmp_sim[j]=chooser_evaluator->multiply_plain((*sim_input_ptr)[j],encodeFractionalChooser(weights[j+(i*in_dim)]));

// 			}
// 		cout<<"b "<<i<<endl;
// 		int bias=biases[i]*approx;
// 		if(bias==0)
// 			tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],31,1);
// 		else
// 			tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],chooser_encoder->encode(bias));
// 		//tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],encodeFractionalChooser(biases[i]));
// 		sim_out[i]=chooser_evaluator->add_many(tmp_sim);
// 	 	}


// 	if(sim_input_ptr!=&sim_input)
// 		delete sim_input_ptr;
	
// 	cout<<" ended fc"<<flush<<endl;
	
// 	return sim_out;
// }

// ChooserPoly FullyConnectedLayer::fullyConnectedSimulator(ChooserPoly sim_input, int in_dim){
// 	cout<<"fc"<<flush;
// 	cout<<"in_dim "<<in_dim<<endl;
// 	vector<ChooserPoly> tmp_sim(in_dim);


// 	sim_input=chooser_evaluator->multiply_plain(sim_input,10,1);

// 	for(int i=0;i<in_dim;i++){
// 		tmp_sim[i]=ChooserPoly(sim_input);
// 	}
	
// 	tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],10,1);

// 	sim_input=chooser_evaluator->add_many(tmp_sim);
	
// 	cout<<" ended fc"<<flush<<endl;
	
// 	return sim_input;
// }


// vector<ChooserPoly> FullyConnectedLayer::fullyConnectedSimulator(vector<ChooserPoly> &sim_input, vector<float> & weights, vector<float> & biases){
// 	cout<<"fc"<<flush;
// 	int approx=1000;
// 	int out_dim=biases.size();
// 	int in_dim=weights.size()/out_dim;
// 	//vector<ChooserPoly> *sim_input_ptr;
// 	vector<ChooserPoly> sim_input_ptr(in_dim);
// 	vector<ChooserPoly> tmp_sim(in_dim);
// 	vector<ChooserPoly> sim_out(1);

// 	if(sim_input.size()<in_dim){
// 		//sim_input_ptr = new vector<ChooserPoly> (in_dim);
// 		int num_copies=in_dim/sim_input.size();
// 		for(int i=0;i<in_dim;i++)
// 			//(*sim_input_ptr)[i]=ChooserPoly(sim_input[int(i/num_copies)]);
// 			sim_input_ptr[i]=ChooserPoly(sim_input[int(i/num_copies)]);

// 	}
// 	//else sim_input_ptr=&sim_input;

// 	for(int i=0;i<1;i++){
// 		for(int j=0; j<in_dim;j++){
// 			cout<<i<<","<<j<<flush<<endl;
// 			int weight=weights[j+(i*in_dim)]*approx;
// 			if(weight==0)
// 				//tmp_sim[j]=chooser_evaluator->multiply_plain((*sim_input_ptr)[j],96,1);
// 				tmp_sim[j]=chooser_evaluator->multiply_plain(sim_input_ptr[j],96,1);
// 			else
// 				//tmp_sim[j]=chooser_evaluator->multiply_plain((*sim_input_ptr)[j],chooser_encoder->encode(weight));
// 				tmp_sim[j]=chooser_evaluator->multiply_plain(sim_input_ptr[j],chooser_encoder->encode(weight));
// 			}
// 		int bias=biases[i]*approx;
// 		if(bias==0)
// 			tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],96,1);
// 		else
// 			tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],chooser_encoder->encode(bias));
// 		sim_out[i]=chooser_evaluator->add_many(tmp_sim);
// 		}


// 	// if(sim_input_ptr!=&sim_input)
// 	// 	delete sim_input_ptr;
	
// 	cout<<" ended fc"<<flush<<endl;
	
// 	return sim_out;
// }






FullyConnectedLayer::~FullyConnectedLayer(){}