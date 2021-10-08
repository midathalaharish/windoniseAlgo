#include <stdio.h>
#include <math.h>
#include "ByquadFilter.h"

int16_t processlpf(BQF *x,int16_t in) {
    float out = in * x->lpf_a[0] + x->lz1;
    x->lz1 = in * x->lpf_a[1] + x->lz2 - x->lpf_b[1] * out;
    x->lz2 = in * x->lpf_a[2] - x->lpf_b[2] * out;
    return out;
}
int16_t processhpf(BQF *x,int16_t in) {
    float out = in * x->lpf_a[0] + x->hz1;
    x->hz1 = in * x->lpf_a[1] + x->hz2 - x->lpf_b[1] * out;
    x->hz2 = in * x->lpf_a[2] - x->lpf_b[2] * out;
    return out;
}
//void GetBEQ(BQF *x)
//{
//    printf("Value of byquad structure Start ** \n");
//    printf("typeEnum :\n bq_type_lowpass = 0, \n bq_type_highpass = 1,\n bq_type_bandpass =2,\n bq_type_notch= 3,\n bq_type_peak =4,\n bq_type_lowshelf=5,\n bq_type_highshelf=6\n\n");
//    printf("ByQuadtype = %d \n", x->type);
//    printf("ByQuadtype Filter taps a0= %lf;a1=%lfa2=%lf;\n b0=1;b1=%lf b2=%lf\n", x->a0,x->a1,x->a2,x->b1,x->b2);
//    printf("ByQuadtype Filter Sampling Rate = %lf\n", x->Fc);
//    printf("ByQuadtype Filter lobe width = %lf\n", x->Q);
//    printf("ByQuadtype Filter PeakGain = %lf\n", x->peakGain);
//    printf("ByQuadtype delay values z1 = %lf , z2 =%lf \n", x->z1, x->z2);
//    printf("Value of byquad structure End ** \n");
//
//}
//void setType(BQF *x, int type) {
//    x->type = type;
//    calcBiquad(x);
//}
//
//void setQ(BQF *x,double Q) {
//    x->Q = Q;
//    calcBiquad(x);
//}
//
//void setFc(BQF *x,double Fc) {
//    x->Fc = Fc;
//    calcBiquad(x);
//}
//
//void setPeakGain(BQF *x,double peakGainDB) {
//    x->peakGain = peakGainDB;
//    calcBiquad(x);
//}
//    
//void setBiquad(BQF *x,int type, double Fc, double Q, double peakGainDB) {
//    x->type = type;
//    x->Q = Q;
//    x->Fc = Fc;
//    setPeakGain(x,peakGainDB);
//    calcBiquad(x);
//}

//void calcBiquad(BQF *x)
//{
//    double norm;
//    double V = pow(10, fabs(x->peakGain) / 20.0);
//    double K = tan(M_PI * x->Fc);
//    switch (x->type) {
//        case bq_type_lowpass:
//            norm = 1 / (1 + K / x->Q + K * K);
//            x->a0 = K * K * norm;
//            x->a1 = 2 * x->a0;
//            x->a2 = x->a0;
//            x->b1 = 2 * (K * K - 1) * norm;
//            x->b2 = (1 - K / x->Q + K * K) * norm;
//            break;
//
//        case bq_type_highpass:
//            norm = 1 / (1 + K / x->Q + K * K);
//            x->a0 = 1 * norm;
//            x->a1 = -2 * x->a0;
//            x->a2 = x->a0;
//            x->b1 = 2 * (K * K - 1) * norm;
//            x->b2 = (1 - K / x->Q + K * K) * norm;
//            break;
//
//        case bq_type_bandpass:
//            norm = 1 / (1 + K / x->Q + K * K);
//            x->a0 = K / x->Q * norm;
//            x->a1 = 0;
//            x->a2 = -(x->a0);
//            x->b1 = 2 * (K * K - 1) * norm;
//            x->b2 = (1 - K / x->Q + K * K) * norm;
//            break;
//
//        case bq_type_notch:
//            norm = 1 / (1 + K / x->Q + K * K);
//            x->a0 = (1 + K * K) * norm;
//            x->a1 = 2 * (K * K - 1) * norm;
//            x->a2 = x->a0;
//            x->b1 = x->a1;
//            x->b2 = (1 - K / x->Q + K * K) * norm;
//            break;
//
//        case bq_type_peak:
//            if (x->peakGain >= 0) {    // boost
//                norm = 1 / (1 + 1/x->Q * K + K * K);
//                x->a0 = (1 + V/x->Q * K + K * K) * norm;
//                x->a1 = 2 * (K * K - 1) * norm;
//                x->a2 = (1 - V/x->Q * K + K * K) * norm;
//                x->b1 = x->a1;
//                x->b2 = (1 - 1/x->Q * K + K * K) * norm;
//            }
//            else {    // cut
//                norm = 1 / (1 + V/x->Q * K + K * K);
//                x->a0 = (1 + 1/x->Q * K + K * K) * norm;
//                x->a1 = 2 * (K * K - 1) * norm;
//                x->a2 = (1 - 1/x->Q * K + K * K) * norm;
//                x->b1 = x->a1;
//                x->b2 = (1 - V/x->Q * K + K * K) * norm;
//            }
//            break;
//        case bq_type_lowshelf:
//            if (x->peakGain >= 0) {    // boost
//                norm = 1 / (1 + sqrt(2) * K + K * K);
//                x->a0 = (1 + sqrt(2*V) * K + V * K * K) * norm;
//                x->a1 = 2 * (V * K * K - 1) * norm;
//                x->a2 = (1 - sqrt(2*V) * K + V * K * K) * norm;
//                x->b1 = 2 * (K * K - 1) * norm;
//                x->b2 = (1 - sqrt(2) * K + K * K) * norm;
//            }
//            else {    // cut
//                norm = 1 / (1 + sqrt(2*V) * K + V * K * K);
//                x->a0 = (1 + sqrt(2) * K + K * K) * norm;
//                x->a1 = 2 * (K * K - 1) * norm;
//                x->a2 = (1 - sqrt(2) * K + K * K) * norm;
//                x->b1 = 2 * (V * K * K - 1) * norm;
//                x->b2 = (1 - sqrt(2*V) * K + V * K * K) * norm;
//            }
//            break;
//        case bq_type_highshelf:
//            if (x->peakGain >= 0) {    // boost
//                norm = 1 / (1 + sqrt(2) * K + K * K);
//                x->a0 = (V + sqrt(2*V) * K + K * K) * norm;
//                x->a1 = 2 * (K * K - V) * norm;
//                x->a2 = (V - sqrt(2*V) * K + K * K) * norm;
//                x->b1 = 2 * (K * K - 1) * norm;
//                x->b2 = (1 - sqrt(2) * K + K * K) * norm;
//            }
//            else {    // cut
//                norm = 1 / (V + sqrt(2*V) * K + K * K);
//                x->a0 = (1 + sqrt(2) * K + K * K) * norm;
//                x->a1 = 2 * (K * K - 1) * norm;
//                x->a2 = (1 - sqrt(2) * K + K * K) * norm;
//                x->b1 = 2 * (K * K - V) * norm;
//                x->b2 = (V - sqrt(2*V) * K + K * K) * norm;
//            }
//            break;
//}
//
//    return;
//}

