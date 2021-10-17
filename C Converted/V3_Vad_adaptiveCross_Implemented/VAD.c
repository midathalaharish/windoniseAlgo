//
//  VAD.c
//  ShivaKumarChikkala
//
//  Created by Midathala Harish on 11.10.21.
//
#include "VAD.h"
void VAD(int16_t *xsample, float E, float tA, int M, int nhop, float vad_threshold, int16_t vad_D[])
{
    int k=0;
    float p = E;
    for(k=0;k<M;k++)
    {
        int16_t c_sample = *(xsample+k);
        p = (1-tA)*p + tA*(c_sample * c_sample);
        if(k==((M/2)-1))
            E=p;
        vad_D[k] = (E >= vad_threshold)? 1:0;
    }
}
