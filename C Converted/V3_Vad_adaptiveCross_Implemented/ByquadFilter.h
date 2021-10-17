#ifndef Biquad_h
#define Biquad_h

#include <stdio.h>
#include <math.h>
#define FTAP_LEN 3
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
    double lpf_b[FTAP_LEN];
    double lpf_a[FTAP_LEN];
    double hpf_b[FTAP_LEN];
    double hpf_a[FTAP_LEN];
    double Fc, Q, peakGain;
    double lz1, lz2;
    double hz1, hz2;
}BQF;

void setType(BQF *x,int type);
void setQ(BQF *x,double Q);
void setFc(BQF *x,double Fc);
void setPeakGain(BQF *x,double peakGainDB);
void Userfeined_Biquad(BQF *x,int type, double Fc, double Q, double peakGainDB);
int16_t processlpf(BQF *x,int16_t in);
int16_t processhpf(BQF *x,int16_t in);
//void calcBiquad(BQF *x);
void GetBEQ(BQF *x);
void SetBEQ(BQF *x);

#endif // Biquad_h
