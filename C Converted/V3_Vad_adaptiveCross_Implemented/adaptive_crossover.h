//
//  adaptive_crossover.h
//  ShivaKumarChikkala
//
//  Created by Midathala Harish on 11.10.21.
//

#ifndef adaptive_crossover_h
#define adaptive_crossover_h

#include <stdio.h>
void adaptive_crossover(int16_t *xl, int16_t *xr, int CutOff_Freq[], uint32_t Fs, float exp_avg_r, int16_t VAD_signal[], uint16_t Win_len, const double *Win, int *Wn);
#endif /* adaptive_crossover_h */
