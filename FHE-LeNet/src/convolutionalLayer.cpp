#include "convolutionalLayer.h"
#include "layer.h"
#include "seal/seal.h"
#include "sealencryptor.h"
#include <cassert>
#include <thread>
#include <ostream>
#include <fstream>

using namespace std;
using namespace seal;


ConvolutionalLayer::ConvolutionalLayer(string name,int xd,int yd,int zd,int xs,int ys,int xf,int yf,int nf, int th_count,plaintext4D & filters, vector<Plaintext> & biases):
	Layer(name), xd(xd), yd(yd), zd(zd),
	xs(xs), ys(ys),
	xf(xf), yf(yf),
	nf(nf),
    /*Number of threads to use to split the computation*/
    th_count(th_count),
	/* Compute the output dimensions of the convolved ciphertext image.
	The function requires the dimensions of the image, the strides, (the padding) and
	the dimensions of the filters. The output is returned through xo,yo and zo. */
	xo( (xd-xf) / xs + 1 ), yo( (yd-yf) / ys + 1 ), zo(  nf),
	filters(filters),
    biases(biases),
    filters_already_ntt(false){
        if(th_count>nf)
            th_count=nf;
    else if(th_count<=0)
           th_count=1;

}
//Initializes weighs and biases by loading them from file_name, where they are saved
ConvolutionalLayer::ConvolutionalLayer(string name,int xd,int yd,int zd,int xs,int ys,int xf,int yf,int nf, int th_count,istream * infile):
    Layer(name), xd(xd), yd(yd), zd(zd),
    xs(xs), ys(ys),
    xf(xf), yf(yf),
    nf(nf),
    th_count(th_count),
    /* Compute the output dimensions of the convolved ciphertext image.
    The function requires the dimensions of the image, the strides, (the padding) and
    the dimensions of the filters. The output is returned through xo,yo and zo. */
    xo( (xd-xf) / xs + 1 ), yo( (yd-yf) / ys + 1 ), zo(  nf),
    filters_already_ntt(false)
    {
        loadPlaintextParameters(infile);
        if(th_count>nf)
            th_count=nf;
    else if(th_count<=0)
           th_count=1;

}

    //da iterare sui 50 kernel
    ciphertext2D ConvolutionalLayer::convolution3d(ciphertext3D image, plaintext3D kernel, Plaintext bias){
        assert(kernel.size()==zd);
        assert(kernel[0].size()==xf);
        assert(kernel[0][0].size()==yf);
        
        int i,j,kx,ky,z,xlast,ylast,p;
        vector<Ciphertext> pixels(xf*yf*zd);
        ciphertext2D result(xo, vector<Ciphertext>(yo));
        Ciphertext image_copy;

        Layer::computeBoundaries(xd,yd,xs,ys,xf,yf, &xlast, &ylast);

        for(i=0;i<xlast;i+=xs){
            assert(i<image[0].size());
            for(j=0;j<ylast;j+=ys){
                assert(j<image[0][0].size());
                p=0;
                for(z=0;z<zd;z++)
                    for(kx=0;kx<xf;kx++)
                        for(ky=0;ky<yf;ky++){
                                //Temporary copy the input ntt Ciphertext
                                image_copy=Ciphertext(image[z][i+kx][j+ky]);
                                //Perform ntt multiplication, with previously ntt transformed kernel
                                evaluator->multiply_plain_ntt(image_copy,kernel[z][kx][ky]);
                                //Transform back the result to a normal Ciphertext
                                evaluator->transform_from_ntt(image_copy);
                                //Insert it in vector pixels
                                pixels[p]=Ciphertext(image_copy);
                                p++;
            }
            //adding bias
            evaluator->add_plain(pixels[0],bias);
            evaluator->add_many(pixels,result[i/xs][j/ys]);
                        }
        }
        return result;

    }
//Transform input in ntt to speedup the multiply_plain, before to start the computation
void ConvolutionalLayer::transform_input_to_ntt(ciphertext3D &input){
    int from=0,to=0;
    int threads=th_count;
    mutex mtx;

    vector<thread> th_vector;


    auto parallelTransform=[&](ciphertext3D &input,mutex &mtx,int from, int to){
       // vector<Ciphertext> tmp(xd*yd*(to-from));
        //int t=0;

        for(int z=from;z<to;z++)
            for(int x=0;x<xd;x++)
                for(int y=0;y<yd;y++){
                    //tmp[t]=Ciphertext(input[z][x][y]);
                    //evaluator->transform_to_ntt(tmp[t]);
                    //t++;
                    evaluator->transform_to_ntt(input[z][x][y]);
                }

       /* while (!mtx.try_lock());
        t=0;
        for(int z=from;z<to;z++)
            for(int x=0;x<xd;x++)
                for(int y=0;y<yd;y++){                
                    input[z][x][y]=Ciphertext(tmp[t]);
                    t++;
                }
        mtx.unlock();*/

    };
    
    if(threads>zd)
        threads=zd;

    int thread_depth=zd/threads;

    for (int i = 0; i < threads; i++){
        from=to;
        if(i<threads-1)
            to+=thread_depth;
        else
            to+=thread_depth + (zd%threads);

        th_vector.emplace_back(parallelTransform, ref(input),ref(mtx),from,to);
    }

    for (int i = 0; i < threads; i++){
        th_vector[i].join();
    }


}

//Transform kernel to ntt to speedup the multiply_plain in forward phase
void ConvolutionalLayer::transform_kernel_to_ntt(int kernel_index){
        for(int z=0;z<zd;z++)
            for(int x=0;x<xf;x++)
                for(int y=0;y<yf;y++)
                    evaluator->transform_to_ntt(filters[kernel_index][z][x][y],MemoryPoolHandle::Global());
}

//Forward with threads
ciphertext3D ConvolutionalLayer::forward (ciphertext3D input){
    int from=0,to=0,thread_kernels=0;
    ciphertext3D convolved(zo,ciphertext2D(xo,vector<Ciphertext>(yo)));
    vector<thread> th_vector;

    auto parallelForward=[&](ciphertext3D &input,ciphertext3D &convolved,int from, int to){
        for(int k=from;k<to;k++){
            if(!filters_already_ntt)
                transform_kernel_to_ntt(k);
            convolved[k]=convolution3d(input,filters[k],biases[k]);
        }

    };
    //Transform input in ntt to speedup the multiply_plain, before to start the computation
    transform_input_to_ntt(input);


    //Each thread will compute the forward on a number of filters
    thread_kernels=nf/th_count;
    
    
    for (int i = 0; i < th_count; i++){
        from=to;
        if(i<th_count-1)
            to+=thread_kernels;
        else
            to+=thread_kernels + (nf%th_count);

        th_vector.emplace_back(parallelForward, ref(input),ref(convolved),from,to);

    }
    for (size_t i = 0; i < th_vector.size(); i++)
    {
        th_vector[i].join();
    }
    
    filters_already_ntt=true;
    return convolved;
}

//Forward without threads
/*
ciphertext3D ConvolutionalLayer::forward (ciphertext3D input){
    Plaintext bias;
	plaintext3D kernel(zd, plaintext2D(xf,vector<Plaintext>(yf)));
    ciphertext3D convolved(zo,ciphertext2D(xo,vector<Ciphertext>(yo)));
	for(int k=0;k<nf;k++){
		kernel=getKernel(k);
        bias=getBias(k);
        convolved[k]=convolution3d(input,kernel,bias);
	}
    return convolved;
}*/

void ConvolutionalLayer::savePlaintextParameters(ostream * outfile){
    int i,j,z,n;

    for(n=0;n<nf;n++){
        for(z=0;z<zd;z++){
            for(i=0;i<xf;i++){
                for(j=0;j<yf;j++){
                    filters[n][z][i][j].save(*outfile);
                    outfile->flush();
                }
            }
        }
        biases[n].save(*outfile);
        outfile->flush();
    }


}
//Load and transform weights in ntt
void ConvolutionalLayer::loadPlaintextParameters(istream * infile){
    int n,z,i,j;
    plaintext4D encoded_weights(nf, plaintext3D(zd, plaintext2D (xf, vector<Plaintext> (yf ) )));
    vector<Plaintext> encoded_biases(nf);

    for(n=0;n<nf;n++){
        for(z=0;z<zd;z++){
            for(i=0;i<xf;i++){
                for(j=0;j<yf;j++){
                    encoded_weights[n][z][i][j].load(*infile);
                }
            }
        }
        encoded_biases[n].load(*infile);
    }
    filters=encoded_weights;
    biases=encoded_biases;
}

plaintext3D ConvolutionalLayer::getKernel(int kernel_index){
    return filters[kernel_index];

}

Plaintext ConvolutionalLayer::getBias(int bias_index){
   return biases[bias_index]; 
 }

void ConvolutionalLayer::printLayerStructure(){
    cerr<<"Convolutional "<<name<<" : input ("<<zd<<","<<xd<<","<<yd<<"); kernel("<<nf<<","<<xf<<","<<yf<<"); stride("<<xs<<","<<ys<<"); output("<<
    zo<<","<<xo<<","<<yo<<") "<<"run with "<<th_count<<" threads"<<endl;
}


// vector<ChooserPoly> ConvolutionalLayer::convolutionalSimulator(vector<ChooserPoly> & sim_input, int xf,int yf, int nf, vector<float> & weights, vector<float> & biases){
//     cout<<"conv"<<flush;
//     int zf=sim_input.size(), kernel_size=xf*yf*zf, approx=1000;
//     vector<ChooserPoly> tmp_sim(kernel_size);
//     vector<ChooserPoly> out_sim(nf);


//     for(int n=0;n<nf;n++){
//         int span=kernel_size*n;
//             for(int i=0;i<xf*yf*zf;i++){
//                     cout<<n<<","<<i<<","<<weights[i+span]*approx<<flush<<endl;
//                     int weight=weights[i+span]*approx;
//                     if(weight==0)
//                         tmp_sim[i]=chooser_evaluator->multiply_plain(sim_input[int(i/(xf*yf))],31,1);
//                     else    
//                         tmp_sim[i]=chooser_evaluator->multiply_plain(sim_input[int(i/(xf*yf))],chooser_encoder->encode(weight));
//                     //tmp_sim[i]=chooser_evaluator->multiply_plain(sim_input[int(i/(xf*yf))],encodeFractionalChooser(weights[i+span]));  
//                 }
//             int bias=biases[n]*approx;
//             if(bias==0)
//                 tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],31,1);
//             else
//                 tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],chooser_encoder->encode(bias));
//             //tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0],encodeFractionalChooser(biases[n]));
//             out_sim[n]=chooser_evaluator->add_many(tmp_sim);
//     }
//     cout<<"ended conv"<<flush;
//     return out_sim;

// }

// ChooserPoly ConvolutionalLayer::convolutionalSimulator(ChooserPoly sim_input, int xf,int yf,int zf){
//     cout<<"conv"<<flush;
//     int kernel_size=xf*yf*zf;
//     vector<ChooserPoly> tmp_sim(kernel_size);
    

//     sim_input=chooser_evaluator->multiply_plain(sim_input,10,1);
    
//     for(int i=0;i<kernel_size;i++){
//         tmp_sim[i]=ChooserPoly(sim_input);

//         }
//     tmp_sim[0]=chooser_evaluator->add_plain(tmp_sim[0], 10,1);

//     sim_input=chooser_evaluator->add_many(tmp_sim);
    
//     cout<<"ended conv"<<flush;
//     return sim_input;

// }



ConvolutionalLayer::~ConvolutionalLayer(){}
