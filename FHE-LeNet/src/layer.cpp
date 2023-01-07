#include "layer.h"


#include <string>

using namespace std;
Layer::Layer(){}
Layer::Layer(string layer_name):name(layer_name){}
Layer::~Layer(){}


void Layer::computeBoundaries(int xd, int yd, int xs, int ys, int xf, int yf, int* xl, int* yl){
    /* The boundary is determined by the maximum value between the filter size
    and the stride */
    if (xf > xs) {
        *xl = xd-xf+1;
    } else {
        *xl = xd-xs+1;
    }
    if (yf > ys) {
        *yl = yd-yf+1;
    } else {
        *yl = yd-ys+1;
    }
    return;

}