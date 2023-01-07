#ifndef UTILS
#define UTILS 

#include <vector>
#include <string>
using namespace std;
//normalizes a dataset (a vector of images) as img'=(img-mean)/std
vector<vector<float> > normalize(vector<vector<float> > dataset, float mean, float stdv);
vector<vector<float> > loadAndNormalizeMNISTestSet(string dataset_path);
vector<unsigned char> loadMNISTestLabels(string dataset_path);
vector<unsigned char> loadMNISTPlainModelPredictions(string file_path);


#endif

