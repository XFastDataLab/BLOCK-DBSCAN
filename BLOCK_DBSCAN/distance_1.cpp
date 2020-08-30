#include "distance_1.h"
#include <iostream>
using namespace std;
float distance_1(float* p1, float* p2)
{
//«–±»—©∑Úæ‡¿Î
    float sum = 0.0;
    for(int i = 0; i < DIME_NUM-1;i++)
    {
       float d =fabsf(p1[i] -p2[i]);
       if(sum<d)
            sum=d;
    }
    return sum;
//≈∑ Ωæ‡¿Î
//    float sum = 0.0;
//    for(int i = 0; i < DIME_NUM-1;i++)
//    {
//        float d1 =p1[i] -p2[i];
//        d1 *= d1;
//        sum = sum+d1;
//    }
//  return sqrt(sum);

}

