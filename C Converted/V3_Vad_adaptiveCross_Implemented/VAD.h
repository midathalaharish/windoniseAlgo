//
//  VAD.h
//  ShivaKumarChikkala
//
//  Created by Midathala Harish on 11.10.21.
//

#ifndef VAD_h
#define VAD_h

#include <stdio.h>
void VAD(int16_t *x_vpu_filtered, float E, float tA, int Win_len, int32_t nhop, float vad_threshold, int16_t vad_D[]);
#endif /* VAD_h */
