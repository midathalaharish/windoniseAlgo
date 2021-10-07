#ifndef Biquad_h
#define Biquad_h

#include <stdio.h>
#include <math.h>

enum {
    bq_type_lowpass = 0,
    bq_type_highpass,
    bq_type_bandpass,
    bq_type_notch,
    bq_type_peak,
    bq_type_lowshelf,
    bq_type_highshelf
};

typedef struct ByQuadFilterStruct{
    int type;
    double a0, a1, a2, b1, b2;
    double Fc, Q, peakGain;
    double z1, z2;
}BQF;

void setType(BQF *x,int type);
void setQ(BQF *x,double Q);
void setFc(BQF *x,double Fc);
void setPeakGain(BQF *x,double peakGainDB);
void Userfeined_Biquad(BQF *x,int type, double Fc, double Q, double peakGainDB);
float process(BQF *x,float in);
void calcBiquad(BQF *x);
void GetBEQ(BQF *x);
void SetBEQ(BQF *x);



#endif // Biquad_h
