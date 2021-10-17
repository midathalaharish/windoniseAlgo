//
//  adaptive_crossover.c
//  ShivaKumarChikkala
//
//  Created by Midathala Harish on 11.10.21.
//

#include "adaptive_crossover.h"
#include "ByquadFilter.h"
#include <math.h>
double FilterBank[][3] ={{2.7214e-05,5.4428e-05,2.7214e-05},{1.0000,-1.9852,0.9853} /* Bl1, Al1*/,
                               {0.9926,-1.9852,0.9926},{1.0000,-1.9852,0.9853} /* Bh1, Ah1*/,
                               {4.1295e-04,8.2590e-04,4.1295e-04},{1.0000,-1.9417,0.9434} /* Bl2, Al2*/,
                               {0.9713,-1.9425,0.9713},{1.0000,-1.9417,0.9434} /* Bh2, Ah2*/,
                               {1.6057e-03,3.2114e-03,1.6057e-03},{1.0000,-1.8835,0.8899} /* Bl3, Al3*/,
                               {0.9434,-1.8867,0.9434},{1.0000,-1.8835,0.8899} /* Bh3, Ah3*/,
                               {5.9885e-03,1.1977e-02,5.9885e-03},{1.0000,-1.7695,0.7934} /* Bl4, Al4*/,
                               {0.8907,-1.7814,0.8907},{1.0000,-1.7695,0.7934} /* Bh4, Ah4*/,
                               {0.021621,0.043241,0.021621},{1.0000,-1.5431,0.6296} /* Bl5, Al5*/,
                               {0.7932,-1.5864,0.7932},{1.0000,-1.5431,0.6296} /* Bh5, Ah5*/,
                               {0.072231,0.144462,0.072231},{1.0000,-1.1092,0.3982} /* Bl6, Al6*/,
                               {0.6268,-1.2537,0.6268},{1.0000,-1.1092,0.3982} /* Bh6, Ah6*/,
};
void adaptive_crossover(int16_t *xframe_l, int16_t *xframe_r, int CutOff_Freq[], uint32_t Fs, float exp_avg_r, int16_t vad_signal[], uint16_t M, const double *Win, int *Wn)
{
    uint16_t k =0;
    BQF BPFLpfHpf;
    int startfilter =0;
    float q=0;
    float alpha_exp = 1/(Fs*0.01);
    for(k=0;k<M;k++)
    {
        if(vad_signal[k] == 0)
        {
            int16_t xframe_r_val = *(xframe_r+k);
            q  =  (1-alpha_exp)*q + alpha_exp * (xframe_r_val*xframe_r_val);
        }
    }
    int32_t noise_est_var_dB  =   (int32_t )(10*log10(q + 1.0e-12));
    if (noise_est_var_dB  <= -60)
    {
        *Wn = CutOff_Freq[0];
    }
    else if(noise_est_var_dB  <= -50)
    {
        *Wn = CutOff_Freq[1];
        startfilter =4;
    }
    else if(noise_est_var_dB  <= -40)
    {
        *Wn = CutOff_Freq[2];
        startfilter =8;
    }
    else if(noise_est_var_dB  <= -30)
    {
        *Wn = CutOff_Freq[3];
        startfilter =12;
    }
    else if(noise_est_var_dB  <= -20)
    {
        *Wn = CutOff_Freq[4];
        startfilter =16;
    }
    else
    {
        *Wn = CutOff_Freq[5];
        startfilter =20;
    }
    BPFLpfHpf.lz1 = BPFLpfHpf.lz2 = 0.0;
    BPFLpfHpf.hz1 = BPFLpfHpf.hz2 = 0.0;
    for(k=0;k<FTAP_LEN;k++)
    {
        BPFLpfHpf.lpf_b[k]= FilterBank[startfilter][k];
        BPFLpfHpf.lpf_a[k]= FilterBank[startfilter+1][k];
    
        BPFLpfHpf.hpf_b[k]= FilterBank[startfilter+2][k];
        BPFLpfHpf.hpf_a[k]= FilterBank[startfilter+3][k];
    }
    
    for(k=0;k<M;k++)
    {
        *(xframe_l + k)  =   (*(Win+k) ) *  (*(xframe_l + k));
        *(xframe_r + k)  =   (*(Win+k) ) *  (*(xframe_r + k));
        
        *(xframe_l + k) = processlpf(&BPFLpfHpf, *(xframe_l+k));
        *(xframe_r + k) = processhpf(&BPFLpfHpf, *(xframe_r+k));
        
        *(xframe_l + k) = *(xframe_l + k) + *(xframe_r + k); //X_com
    }

}
