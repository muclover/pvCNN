#include "sealencryptor.h"
#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"
#include <fstream>
#include <string>
#include <ostream>
#include <cassert>

using namespace seal;
using namespace std;

EncryptionParameters * parms;
SEALContext * context;
KeyGenerator * keygen;
Encryptor * encryptor;
Decryptor * decryptor;
Evaluator * evaluator;
IntegerEncoder * intencoder;
FractionalEncoder * fraencoder;
EvaluationKeys * ev_keys16;

ChooserEncoder * chooser_encoder;
ChooserEvaluator * chooser_evaluator;

void setParameters(int poly_modulus, uint64_t plain_modulus){

  //4096,1<<20
	parms = new EncryptionParameters();
	parms->set_poly_modulus("1x^"+to_string(poly_modulus)+" + 1");
  parms->set_coeff_modulus(coeff_modulus_128(poly_modulus));
  parms->set_plain_modulus(plain_modulus);


    context = new SEALContext(*parms);

    keygen= new KeyGenerator(*context);
    auto public_key = keygen->public_key();
    auto secret_key= keygen->secret_key();
    encryptor = new Encryptor(*context, public_key);
    decryptor = new Decryptor(*context,secret_key);
    evaluator = new Evaluator(*context);
    //intencoder = new IntegerEncoder(context->plain_modulus(),3);
    //64,32
    fraencoder =  new FractionalEncoder(context->plain_modulus(), context->poly_modulus(), 64, 32, 3);
   	ev_keys16 = new EvaluationKeys();
    keygen->generate_evaluation_keys(16, *ev_keys16);
   
}

void setAndSaveParameters(string public_key_path,string secret_key_path,string evaluation_key_path, int poly_modulus, uint64_t plain_modulus){
    ofstream pubfile(public_key_path, ofstream::binary);
    ofstream secfile(secret_key_path, ofstream::binary);
    ofstream evalfile(evaluation_key_path, ofstream::binary);

    setParameters(poly_modulus,plain_modulus);

    keygen->public_key().save(pubfile);
    keygen->secret_key().save(secfile);
    ev_keys16->save(evalfile);
    
    pubfile.close();
    secfile.close();
    evalfile.close();


}


void initFromKeys(string public_key_path,string secret_key_path,string evaluation_key_path,int poly_modulus, uint64_t plain_modulus){
	ifstream pubfile(public_key_path, ifstream::binary);
  ifstream secfile(secret_key_path, ifstream::binary);
  ifstream evalfile(evaluation_key_path, ifstream::binary);

	parms = new EncryptionParameters();
	parms->set_poly_modulus("1x^"+to_string(poly_modulus)+" + 1");
  parms->set_coeff_modulus(coeff_modulus_128(poly_modulus));
  parms->set_plain_modulus(plain_modulus);
  context = new SEALContext(*parms);


  PublicKey public_key;
  public_key.load(pubfile);
  SecretKey secret_key;
  secret_key.load(secfile);

  keygen=new KeyGenerator(*context,secret_key,public_key);

  encryptor = new Encryptor(*context, public_key);
  decryptor = new Decryptor(*context,secret_key);
  evaluator = new Evaluator(*context);
  //intencoder = new IntegerEncoder(context->plain_modulus(),3);
  fraencoder =  new FractionalEncoder(context->plain_modulus(), context->poly_modulus(), 64, 32, 3);


  ev_keys16 = new EvaluationKeys();
  ev_keys16->load(evalfile);

  pubfile.close();
  secfile.close();
  evalfile.close();
  
   
}


void delParameters(){
   delete ev_keys16;
   delete fraencoder;
   delete intencoder;
   delete evaluator;
   delete decryptor;
   delete encryptor;
   delete keygen;
   delete context;
   delete parms;
}
//Precondition: setParameters() or initFromKeys() must be called before 
//Encrypt a normalized image(that is a vector of float) using Fractional encoder in a 3d ciphertext of dim zd,xd,yd
ciphertext3D encryptImage(vector<float> image, int zd, int xd, int yd){

  ciphertext3D encrypted_image(zd,ciphertext2D(xd,vector<Ciphertext>(yd)));
  
    for(int z=0;z<zd;z++)
      for(int i=0;i<xd;i++)
        for(int j=0;j<yd;j++){
            cout<<"Image to be encrypted:"<<image[i*xd+j]<<"!"<<endl;
            encryptor->encrypt(fraencoder->encode(image[i*xd+j]),encrypted_image[z][i][j]);          
            }

  cout << "Noise budget in fresh encrypted pixel: ";
  //cout<< decryptor->invariant_noise_budget(encrypted_image[0][0][0]) << " bits"<<endl;

  return encrypted_image;

}

ciphertext3D encryptImage(floatCube image){

  int zd=image.size(), xd=image[0].size(), yd=image[0][0].size();
  ciphertext3D encrypted_image(zd,ciphertext2D(xd,vector<Ciphertext>(yd)));

  for(int z=0;z<zd;z++)
    for(int i=0;i<xd;i++)
      for(int j=0;j<yd;j++){
          encryptor->encrypt(fraencoder->encode(image[z][i][j]),encrypted_image[z][i][j]);          
          }

  return encrypted_image;

}

// ciphertext3D deepCopyImage(ciphertext3D image){
//     int zd=image.size(), xd=image[0].size(), yd=image[0][0].size();
//     ciphertext3D encrypted_image(zd,ciphertext2D(xd,vector<Ciphertext>(yd)));

//     for(int z=0;z<zd;z++)
//       for(int i=0;i<xd;i++)
//         for(int j=0;j<yd;j++){
//           encrypted_image[z][i][j]=Ciphertext(image[z][i][j]);         
//           }

//   return encrypted_image;
// }


//Precondition: setParameters() or initFromKeys() must be called before
//Encrypt a normalized image(that is a vector of float) using Fractional encoder in a 3d ciphertext of dim zd,xd,yd and save it in file file_name
ciphertext3D encryptAndSaveImage(vector<float> image, int zd, int xd, int yd, string file_name){
  ofstream outfile(file_name, ofstream::binary);
  ciphertext3D encrypted_image(zd,ciphertext2D(xd,vector<Ciphertext>(yd)));
  
    for(int z=0;z<zd;z++)
      for(int i=0;i<xd;i++)
        for(int j=0;j<yd;j++){
            encryptor->encrypt(fraencoder->encode(image[i*xd+j]),encrypted_image[z][i][j]);
            encrypted_image[z][i][j].save(outfile);       
            }

  outfile.close();
  cout << "Noise budget in fresh encrypted pixel: ";
  cout<< decryptor->invariant_noise_budget(encrypted_image[0][0][0]) << " bits"<<endl;

  return encrypted_image;
}




//Precondition: initFromKeys() must be called before and encryptAndSaveImage() with that parameters must have been called before
//Load a normalized image saved it in file file_name in a 3d ciphertext of dim zd,xd,yd and 
ciphertext3D loadEncryptedImage(int zd, int xd, int yd, string file_name){
    ifstream imagefile(file_name, ifstream::binary);
    ciphertext3D encrypted_image(zd,ciphertext2D(xd,vector<Ciphertext>(yd)));

    for(int z=0;z<zd;z++)
      for(int i=0;i<xd;i++)
        for(int j=0;j<yd;j++){
            encrypted_image[z][i][j].load(imagefile);       
            }

    imagefile.close();
    return encrypted_image;
}
//Decrypt a fractional encrypted 3d image
floatCube decryptImage(ciphertext3D encrypted_image){
  		
    int zd=encrypted_image.size();
    int xd=encrypted_image[0].size();
    int yd=encrypted_image[0][0].size();
    //cerr<<zd<<","<<xd<<","<<yd<<endl;

    Plaintext tmp;
    floatCube image(zd, vector<vector<float> > (xd,vector<float>(yd)));

    for(int z=0;z<zd;z++){
        //cerr<<z<<endl;
        for(int i=0;i<xd;i++){
            for(int j=0;j<yd;j++){
                    decryptor->decrypt(encrypted_image[z][i][j], tmp);
                    image[z][i][j]=fraencoder->decode(tmp);
                    //cout<<image[z][i][j]<<",";
        }
    //cout<<endl;
    }
    //cout<<endl;
    //cout<<endl;
  }
  cout<<"Decrypt Image Success!"<<endl;
  return image;
}

void setChooserParameters(uint64_t base){
  chooser_encoder = new ChooserEncoder(base);
  //fraencoder =  new FractionalEncoder(SmallModulus(base), BigPoly("1x^1024 + 1"), 64, 32, 3);
  chooser_evaluator = new ChooserEvaluator;

}

void delChooserParameters(){
  delete chooser_evaluator;
  //delete fraencoder;
  delete chooser_encoder;
}

ChooserPoly encodeFractionalChooser(float value){
  ChooserPoly c_poly;
  Plaintext value_poly;
  value_poly=fraencoder->encode(value);
  c_poly.reset();
  c_poly.max_coeff_count() = value_poly.significant_coeff_count();
  c_poly.max_abs_value() = util::poly_infty_norm_coeffmod(
  value_poly.data(), value_poly.coeff_count(), fraencoder->plain_modulus());
  return c_poly;
}


/*
Helper function: Prints the parameters in a SEALContext.
*/
void print_parameters(const SEALContext &context)
{
    cout << "/ Encryption parameters:" << endl;
    cout << "| poly_modulus: " << context.poly_modulus().to_string() << endl;

    /*
    Print the size of the true (product) coefficient modulus
    */
    cout << "| coeff_modulus size: " 
        << context.total_coeff_modulus().significant_bit_count() << " bits" << endl;

    //Number of qi in coeff_count 
    int coeff_mod_count = static_cast<int>(context.coeff_modulus().size());

    cout << "| coeff_modulus qi: ";
    for (int i = 0; i < coeff_mod_count; i++)
      cout<<"coeff "<<context.coeff_modulus()[i].value();
      cout<<endl;


    cout << "| plain_modulus: " << context.plain_modulus().value() << endl;
    cout << "\\ noise_standard_deviation: " << context.noise_standard_deviation() << endl;
    cout << endl;
}